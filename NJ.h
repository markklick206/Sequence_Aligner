#include "NW_MultiAlign.h"
//#define PRINTDISTMAT
//#define PRINTQMAT
//#define PRINTNEWDISTMAT

double** CreateDistanceMatrix(std::vector<MultiSequence*> *SeqSet) {
	NWMultiAlign* Align = new NWMultiAlign();

	int x = SeqSet->size();

	double** DistMat = new double*[x];
	for (int i = 0; i < x; i++)
		DistMat[i] = new double[x];

	for (int i = 0; i < x; i++) {
		for (int j = i; j < x; j++) {
			Align->SetMultiSequence((*SeqSet)[i], 1);
			Align->SetMultiSequence((*SeqSet)[j], 2);
			Align->AlignMultiSequences();
			DistMat[i][j] = Align->LevenshteinDistance();
			DistMat[j][i] = DistMat[i][j];
		}
	}
	return DistMat;
}

double Q(double** DistMat, int u, int v, int n) {
	double q;
	double sumik = 0, sumjk = 0;

	for (int k = 0; k < n; k++)
		sumik += DistMat[u][k];

	for (int k = 0; k < n; k++)
		sumjk += DistMat[v][k];

	q = (n - 2) *DistMat[u][v] - sumik - sumjk;

	return q;
}

double d(double** DistMat, int k, int f, int g) {
	return 0.5 * (DistMat[f][k] + DistMat[g][k] - DistMat[f][g]);
}

void NeighborJoin(std::vector<MultiSequence*> *SeqSet) {
	double** DistMat = CreateDistanceMatrix(SeqSet);

	while (SeqSet->size() > 2) {
		int x = SeqSet->size();

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

		double** QMat = new double*[x];
		for (int i = 0; i < x; i++)
			QMat[i] = new double[x];


		double minQ = (int)INFINITY;
		int SeqI, SeqJ;
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < x; j++) {
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

		NWMultiAlign* Aligner = new NWMultiAlign();
		Aligner->SetMultiSequence((*SeqSet)[SeqI], 1);
		Aligner->SetMultiSequence((*SeqSet)[SeqJ], 2);

		Aligner->AlignMultiSequences();

		SeqSet->push_back(Aligner->GetAlignedMultiSequence());

		if (SeqI > SeqJ) {
			SeqSet->erase(SeqSet->begin() + SeqI);
			SeqSet->erase(SeqSet->begin() + SeqJ);
		}
		else {
			SeqSet->erase(SeqSet->begin() + SeqJ);
			SeqSet->erase(SeqSet->begin() + SeqI);
		}

		x = SeqSet->size();
		double** newDistMat = new double*[x];
		for (int i = 0; i < x; i++)
			newDistMat[i] = new double[x];

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
		for (int i = 0; i < x; i++) {
			if (SeqI == oldI || SeqJ == oldI)
				oldI++;
			if (SeqI == oldI || SeqJ == oldI)
				oldI++;
			newDistMat[i][x - 1] = d(DistMat, oldI, SeqI, SeqJ);
			oldI++;
		}

		oldJ = 0;
		for (int i = 0; i < x; i++) {
			if (SeqI == oldJ || SeqJ == oldJ)
				oldJ++;
			if (SeqI == oldJ || SeqJ == oldJ)
				oldJ++;
			newDistMat[x - 1][i] = d(DistMat, oldJ, SeqI, SeqJ);
			oldJ++;
		}
		newDistMat[x - 1][x - 1] = 0;

		DistMat = newDistMat;

		//std::cout << SeqSet->size() << std::endl;
		//for (unsigned int i = 0; i < SeqSet->size(); i++) {
		//	std::string str;
		//	str = "ms";
		//	str.push_back(65 + i + rand()%23);
		//	str.append(".txt");
		//	(*SeqSet)[i]->WriteMultiSequenceToFile(str);
		//}
	}

	NWMultiAlign* Aligner = new NWMultiAlign();
	Aligner->SetMultiSequence((*SeqSet)[0], 1);
	Aligner->SetMultiSequence((*SeqSet)[1], 2);

	Aligner->AlignMultiSequences();
	//Aligner->WriteAlignedMultiSequenceToFile("AAAA.txt");

	SeqSet->erase(SeqSet->begin() + 1);
	SeqSet->erase(SeqSet->begin());

	SeqSet->push_back(Aligner->GetAlignedMultiSequence());

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
}





