#ifndef SAVECLASS_H
#define SAVECLASS_H

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <string>
#include <vector>
#include <stdlib.h>

#include "../../../baseTools/tools/constants.h"
#include "atomClass.h"
#include "moleculeClass.h"
#include "pdfClass.h"
#include "diffractionClass.h"

using namespace std;


class SAVEclass {

   public:
	SAVEclass();
	SAVEclass(string filename, string treename);
	~SAVEclass();

	TFile* file;
	TTree* tree;

        void newFile(string filename, string treename);
	void addInfo(DIFFRACTIONclass* diff_in);
	void addInfo(vector<string> names, vector<double*> ptrs);
	void fillEvent();
	void writeFile();

   private:

	DIFFRACTIONclass* diff;

};

#endif
