//
//  processvcfstats.h
//  VCF_Comparisons
//
//  Created by Pelin Sahlen on 17/09/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//
struct stats{
    double x[8][2];
    
 //   double nof_snps[2];
  //  double nof_dbsnps[2];
   // double nof_novelsnps[2];
    //double trtvratio[2];
   // double homhetratio[2];
    double discordantgens;
    double discordantdbsnpgens;
   // double nofcommonsnps[2];
   // double nofcommondbsnps[2];
   // double nofcommonnovelsnps[2];
    double nofvalsnps;
    double valsnpsfoundinboth;
    double valsnpsfoundinfirst;
    double valsnpsfoundinsecond;
};

vector < stats > statsvector;

void readstatfile(string filename){
    statsvector.push_back(stats());
    ifstream infile(filename.c_str());
    string temp;
    getline(infile,temp);
    getline(infile,temp);
    infile >> statsvector.back().x[0][0] >> statsvector.back().x[1][0] >> statsvector.back().x[2][0];
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    infile >> temp >> statsvector.back().x[3][0] >> temp >> statsvector.back().x[4][0];
    
    getline(infile,temp);
    getline(infile,temp);
    infile >> statsvector.back().x[0][1] >> statsvector.back().x[1][1] >> statsvector.back().x[2][1];
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    getline(infile,temp);
    infile >> temp >> statsvector.back().x[3][1] >> temp >> statsvector.back().x[4][1];
    
    getline(infile,temp);
    infile >> temp >> statsvector.back().nofvalsnps >> temp >> temp >> statsvector.back().valsnpsfoundinfirst
           >> temp >> temp >> statsvector.back().valsnpsfoundinsecond
           >> temp >> temp >> statsvector.back().valsnpsfoundinboth;
    
    getline(infile,temp);
    infile >> temp >> temp >> temp >> statsvector.back().discordantgens >> temp >> temp >> temp >> statsvector.back().discordantdbsnpgens;
    getline(infile,temp);
    infile >> temp >> temp >> temp >> temp >> temp >> statsvector.back().x[5][0] >> statsvector.back().x[5][1];
    infile >> temp >> temp >> temp >> temp >> temp >> statsvector.back().x[6][0] >> statsvector.back().x[6][1];
    infile >> temp >> temp >> temp >> temp >> temp >> statsvector.back().x[7][0] >> statsvector.back().x[7][1];
    
    infile.close();
    
}
void readvcfstats(void){
    
    vector< string > filenames;
    string fname1, fname2,buffer;
    
    ofstream outfile("stats_averaged.txt");
    
    ifstream filelist("filelist.txt");

    filelist >> fname1 >> fname2;
    int a = fname1.find("/");
    int b = fname2.find("/");
    string comparison1 = fname1.substr(0,a);
    string comparison2 = fname2.substr(0,b);
    cout << fname1 << "  " << fname2 << "  " << comparison1 << "  " << comparison2 << endl;
    
do{
    while (comparison1 == comparison2){
        filenames.push_back(fname1);
        filenames.push_back(fname2);
        fname1 = fname2;
        comparison1 = comparison2;
        filelist >> fname2;
        int b = fname2.find("/");
        string comparison2 = fname2.substr(0,b);
        cout << "h " << fname2 << " " << comparison2 << endl;
    }
    
    for (int i = 0; i < filenames.size(); ++i) {
        readstatfile(filenames[i]);
    }
    
    double ms1[8], ms2[8];
    double ss1[8], ss2[8];

    alglib::real_1d_array d1;
    alglib::real_1d_array d2;
    std::vector<double> y1, y2;
    
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < filenames.size(); ++i) {
            y1.push_back(statsvector[i].x[i][0]);
            y1.push_back(statsvector[i].x[i][1]);
        }
        d1.setcontent(y1.size(), &(y1[0]));
        d2.setcontent(y2.size(), &(y2[0]));
        
        ms1[i] = alglib::samplemean(d1);
        alglib::sampleadev(d1, ss1[i]);
        
        ms2[i] = alglib::samplemean(d2);
        alglib::sampleadev(d2, ss2[i]);

        outfile << i << '\t' << ms1[i] << " ± " << ss1[i] << '\t' << ms2[i] << " ± " << ss2[i] << endl;

    outfile << "Number of discordant snps" << '\t' << statsvector[i].discordantgens << endl
            << "Number of discordant dbsnps" << '\t' << statsvector[i].discordantdbsnpgens << endl
            << "Number of val snps" << '\t' << statsvector[i].nofvalsnps << endl
            << "valsnpsfoundinboth" << '\t' << statsvector[i].valsnpsfoundinboth << endl
            <<  "valsnpsfoundinfirst" << '\t' << statsvector[i].valsnpsfoundinfirst << endl
            <<  "valsnpsfoundinsecond" << '\t' << statsvector[i].valsnpsfoundinsecond << endl;

    }
    statsvector.resize(0);
    filenames.resize(0);
    fname1 = fname2;
    comparison1 = comparison2;
    filelist >> fname2;
    int b = fname2.find("/");
    string comparison2 = fname2.substr(0,b);
    
}while(fname2 != "END");
}
