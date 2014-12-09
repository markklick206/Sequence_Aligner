#include "NW_MultiAlign.h"
//#define PRINTDISTMAT
//#define PRINTQMAT
//#define PRINTNEWDISTMAT

VVD CreateDistanceMatrix(std::vector<MultiSequence>& SeqSet) {
	int x = static_cast<int>(SeqSet.size());
	VVD DistMat;
	DistMat.resize(x);
	for (int i = 0; i < x; i++)
		DistMat[i].resize(x);

	NWMultiAlign Align;

	std::cout << "Starting to create DistMat..." << std::endl;
	int progress = 0, load = x * (x + 1) / 2;
	for (int i = 0; i < x; i++) {
		Align.SetMultiSequence(SeqSet[i], 1);
		for (int j = i; j < x; j++) {
			//std::cout << "At line " << __LINE__ << " " << i << " " << j << std::endl;
			Align.SetMultiSequence(SeqSet[j], 2);
			Align.AlignMultiSequences();
			DistMat[i][j] = Align.LevenshteinDistance();
			Align.ClearMSOut();
			DistMat[j][i] = DistMat[i][j];
			progress++;
			std::cout << "Working -- " << progress << "/" << load << std::endl;
		}
	}
	return DistMat;
}

double Q(VVD& DistMat, int u, int v, int n) {
	double q;
	double sumik = 0, sumjk = 0;

	for (unsigned int k = 0; k < DistMat.size(); k++)
		sumik += DistMat[u][k];

	for (unsigned int k = 0; k < DistMat.size(); k++)
		sumjk += DistMat[v][k];

	q = (n - 2) * DistMat[u][v] - sumik - sumjk;

	return q;
}

double d(VVD& DistMat, int k, int f, int g) {
	return 0.5 * (DistMat[f][k] + DistMat[g][k] - DistMat[f][g]);
}

void NeighborJoin(std::vector<MultiSequence>& SeqSet) {
	
	VVD DistMat = CreateDistanceMatrix(SeqSet);
	VVD QMat;
	int x = 0;

	while (SeqSet.size() > 2) {
		x = static_cast<int>(SeqSet.size());
        /*
#ifdef PRINTDISTMAT
		std::ofstream output;
		std::stringstream ss;
		LARGE_INTEGER time;
		QueryPerformanceCounter(&time);
		int a = static_cast<int>(time.QuadPart);
		std::string filename;
		ss << "DistMats" << a << ".txt";
		filename = ss.str();
		output.open(filename);
		output << std::endl;
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < x; j++) {
				output << DistMat[i][j] << "\t";
			}
			output << std::endl;
		}
		output << std::endl << std::endl;
		output.close();
#endif
         */

		QMat.resize(x);
		for (int i = 0; i < x; i++)
			QMat[i].resize(x);
		
		
		double minQ = (int)INFINITY;
		int SeqI, SeqJ;
		for (unsigned int i = 0; i < QMat.size(); i++) {
			for (unsigned int j = 0; j < QMat.size(); j++) {
				
				if (i != j)
					QMat[i][j] = Q(DistMat, i, j, x);
				else
					QMat[i][j] = 0;
				
				if (QMat[i][j] < minQ) {
					SeqI = i;
					SeqJ = j;
					minQ = QMat[i][j];
				}
			}
		}
		
        /*
#ifdef PRINTQMAT
		std::ofstream output1;
		std::stringstream ss1;
		time;
		QueryPerformanceCounter(&time);
		a = static_cast<int>(time.QuadPart);
		filename.clear();
		ss1 << "QMats" << a << ".txt";
		filename = ss1.str();
		output1.open(filename);
		output1 << std::endl;
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < x; j++) {
				output1 << QMat[i][j] << "\t";
			}
			output1 << std::endl;
		}
		output1.close();
#endif
         */

		NWMultiAlign* Aligner = new NWMultiAlign();
		Aligner->SetMultiSequence(SeqSet[SeqI], 1);
		Aligner->SetMultiSequence(SeqSet[SeqJ], 2);

		Aligner->AlignMultiSequences();

		SeqSet.push_back(Aligner->GetAlignedMultiSequence());

		Aligner->~NWMultiAlign();
        
		delete Aligner;
		
		if (SeqI > SeqJ) {
			SeqSet.erase(SeqSet.begin() + SeqI);
			SeqSet.erase(SeqSet.begin() + SeqJ);
		}
		else {
			SeqSet.erase(SeqSet.begin() + SeqJ);
			SeqSet.erase(SeqSet.begin() + SeqI);
		}
		
		x = static_cast<int>(SeqSet.size());
        
		VVD newDistMat;
		newDistMat.resize(x);
		for (int i = 0; i < x; i++)
			newDistMat[i].resize(x);

		int oldI = 0, oldJ = 0;
		
		for (int i = 0; i < x - 1; i++) { //JOHN
			if (SeqI == oldI || SeqJ == oldI)
				oldI++;
			if (SeqI == oldI || SeqJ == oldI)
				oldI++;
			oldJ = 0;
			for (int j = 0; j < x - 1; j++) {
				if (SeqI == oldJ || SeqJ == oldJ)
					oldJ++;
				if (SeqI == oldJ || SeqJ == oldJ)
					oldJ++;
				newDistMat[i][j] = DistMat[oldI][oldJ];
				oldJ++;
			}
			oldI++;
		}
		
		oldI = 0;
		for (int i = 0; i < x - 1; i++) {
			if (SeqI == oldI || SeqJ == oldI)
				oldI++;
			if (SeqI == oldI || SeqJ == oldI)
				oldI++;
			newDistMat[i][x - 1] = d(DistMat, oldI, SeqI, SeqJ);
			oldI++;
		}
		
		oldJ = 0;
		for (int i = 0; i < x - 1; i++) {
			if (SeqI == oldJ || SeqJ == oldJ)
				oldJ++;
			if (SeqI == oldJ || SeqJ == oldJ)
				oldJ++;
			newDistMat[x - 1][i] = d(DistMat, oldJ, SeqI, SeqJ);
			oldJ++;
		}
		newDistMat[x - 1][x - 1] = 0;
		
		DistMat = newDistMat;
	}
	
	NWMultiAlign* Aligner = new NWMultiAlign();
	Aligner->SetMultiSequence(SeqSet[0], 1);
	Aligner->SetMultiSequence(SeqSet[1], 2);

	Aligner->AlignMultiSequences();

	SeqSet.erase(SeqSet.begin() + 1);
	SeqSet.erase(SeqSet.begin());

	SeqSet.push_back(Aligner->GetAlignedMultiSequence());

	Aligner->~NWMultiAlign();
    
	delete Aligner;
    /*
#ifdef PRINTNEWDISTMAT
	std::ofstream output1;
	output1.open("newDistMat.txt");
	output1 << std::endl;
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < x; j++) {
			output1 << newDistMat[i][j] << "\t";
		}
		output1 << std::endl;
	};
	output1.close();
#endif
     */
}





