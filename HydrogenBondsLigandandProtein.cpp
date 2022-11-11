#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <bits/stdc++.h>
#include <time.h>
#define pi 3.14159

//To calculate Drug Donor - Drug Donor Hydrogen / Acc Ox - Protein Don/Acceptor HBonds (Explicit Solvent)
// ms
using namespace std;

int nframes,framestart,drugrmsd,nligatoms,nligands;
//frames
string line;
string ATOM,atnum,atname,resname,chain,resid,x,y,z,t,iframe;
//info file
string infofile;
string structure,traj,proteinchain,selection,proteinid,Proteinmain,Proteinside,atindex,rmsd,mdsys,inputfile,inputdata;
double xl,yl,zl;
//Atom index file
string dhatom,drugdindex,drugdhindex;

//vectors
//index
vector <int> drugdonorindex,drugdonorhydindex,proresid,Totalmain,Totalside,mainoratomnum,sideoratomnum,proteinmainind,proteinsideind,mainind,sideind,indx,totalh,countwoangl;
vector <string> dname,res,atoms;
string mainatoms,moind,macind,resn;
double dotp,magddh,magdhp,radangle,angle;
string sideat,sliatom,s_or_atomn,mainat,miatom,m_or_atomn,residue;
//coordinates
vector <double> pcoord,dcoord,hdist,vcddh,vcdhp;
vector <vector<double> > procoord,drugcoord,hbdata;
vector <string> datain;


//getinfo
//
int getinfo(string datafile)
    {
    char data[datafile.length()];
    strcpy(data,datafile.c_str());
    ifstream file;
    file.open(data);
    if (file.is_open())
       {
       while(getline(file,line))
            {
            datain.push_back(line);
            }
       }
    }

//drugdonorindexfile
int drug_index(string indexfile)
   {
   int n=indexfile.length();
   char filearray[n];
   strcpy(filearray,indexfile.c_str());
   ifstream file;
   file.open(filearray);
   if (file.is_open())
      {  
      while(getline(file,line))
          {
          istringstream iss(line);
          iss >> dhatom >> drugdindex >> drugdhindex;
          drugdonorindex.push_back(::atoi(drugdindex.c_str()));
          drugdonorhydindex.push_back(::atoi(drugdhindex.c_str()));
          }
     }
   }

void convert (string instring, ofstream& s)
  {
    int m=instring.length();
    char mc[m];
    strcpy(mc,instring.c_str());
    s.open(mc);
  }


//frame files
int catdcd(string psffile,string dcdfile,int i)
   {
   ostringstream oss;
   oss << "catdcd -o frame.pdb -otype pdb -stype psf -s " << psffile <<" -first " << i << " -last " << i << " -dcd " << dcdfile;
   string comm=oss.str();
   char com[comm.length()];
   strcpy(com, comm.c_str());
   system(com);
   cout << "frame written" << endl;
   }
//

//getcoord
//read frames get drug and protein coord
int drug_getcoord(string filename)    //filename
    {                                                         
    int n=filename.length();                                           
    char filearray[n];
    strcpy(filearray,filename.c_str());
    ifstream file;
    file.open(filearray);
    if (file.is_open())
       {
       while(getline(file,line))
           {
           istringstream iss(line);
           iss>>ATOM>>atnum>>atname>>resname>>resid>>x>>y>>z>>t;           
           if ((resname[0] == 'S') && (resname[1] == 'C') && (resname[2] == 'B'))
              {
              dcoord.push_back(::atof(x.c_str()));
              dcoord.push_back(::atof(y.c_str()));
              dcoord.push_back(::atof(z.c_str()));
              dcoord.push_back(::atoi(atnum.c_str()));
              drugcoord.push_back(dcoord);
              dcoord.clear();
              dname.push_back(atname);
              }
           }//getline
       }// open
    } //drug

int protein_getcoord(string filename,string pchain)
   {
    int n=filename.length();
    char filearray[n];
    strcpy(filearray,filename.c_str());
    ifstream file;
    file.open(filearray);
    if (file.is_open())
       {
       while(getline(file,line))
           {
           istringstream iss(line);
           iss>>ATOM>>atnum>>atname>>resname>>chain>>resid>>x>>y>>z>>t;
           if ( (chain == pchain) && (ATOM == "ATOM"))
              {
              pcoord.push_back(::atof(x.c_str()));
              pcoord.push_back(::atof(y.c_str()));
              pcoord.push_back(::atof(z.c_str()));
              procoord.push_back(pcoord);
              pcoord.clear();
              proresid.push_back(::atoi(resid.c_str()));
              res.push_back(::resname);  //resname
              indx.push_back(::atoi(atnum.c_str()));
	      atoms.push_back(::atname);
              }
           }
       }
     }//protein


//side atoms data
//HN      470     2       N       469     1       SER
int proteinsidemain(string proteinsidemainindexfile)
   {
   int mp=proteinsidemainindexfile.length();
   char filearraymp[mp];
   strcpy(filearraymp,proteinsidemainindexfile.c_str());
   ifstream filemp;
   filemp.open(filearraymp);
   if (filemp.is_open())
      {
      while(getline(filemp,line))
           {
           istringstream iss(line); 
           iss >> sideat >> sliatom >> s_or_atomn >> mainat >> miatom >> m_or_atomn >> residue;
           sideoratomnum.push_back(::atoi(sliatom.c_str()));    //o
           proteinsideind.push_back(::atoi(s_or_atomn.c_str()));
           }
      }
   }

//main atoms data
int proteinmain(string matoms_in)
   {
   int mains=matoms_in.length();
   char datain[mains];
   strcpy(datain,matoms_in.c_str());
   ifstream infile;
   infile.open(datain);
   if (infile.is_open())
      {
      while(getline(infile,line))
           {
           istringstream iss(line);
           iss >> mainatoms >> moind >> macind >> resn;
           mainoratomnum.push_back(::atoi(moind.c_str()));
           proteinmainind.push_back(::atoi(macind.c_str()));
// sort mainat and delete repetitions
//           sort(proteinmainind.begin(), proteinmainind.end());
//           proteinmainind.erase( unique( proteinmainind.begin(), proteinmainind.end() ), proteinmainind.end() );
           }
      }
   }

//Drugmain - Drug - Protein side
double getdist_psidem_druglig(vector <double> vecDH, vector <double> vecPS, int ligid, int dm_id, int dho_id, int p_oid, int p_aid, double x, double y, double z)
    //d p lid mid hoid poid pacid x y z
    {
    double dcal,dist,xcomp,ycomp,zcomp,xc,yc,zc,ccount=0;
    dcal = 0.0;
    dist = 0.0;
    xcomp = vecDH[0] - vecPS[0];
    ycomp = vecDH[1] - vecPS[1];
    zcomp = vecDH[2] - vecPS[2];
    if (abs(xcomp) <= x/2)             xc = abs(xcomp);
    if (abs(ycomp) <= y/2)             yc = abs(ycomp);
    if (abs(zcomp) <= z/2)             zc = abs(zcomp);
    if (abs(xcomp) > (x/2))             xc = x-xcomp;
    if (abs(ycomp) > (y/2))             yc = y-ycomp;
    if (abs(zcomp) > (z/2))             zc = z-zcomp;
    dcal = dcal + (xc * xc) + (yc * yc) + (zc * zc);
    dist = sqrt(dcal);
    if ((dist < 3.00) && (dist != 0.00))
       {
       ccount = ccount + 1;
       hdist.push_back(ligid); //l
       hdist.push_back(dm_id); //dmain id
       hdist.push_back(dho_id);  //dh id
       hdist.push_back(p_oid);   //original protein id   -   as in pbd file
       hdist.push_back(p_aid);   // actual protein id   -   as in vector
       hdist.push_back(dist);  // dist
       hbdata.push_back(hdist);
       hdist.clear();
       }
   return ccount;
   }

/*//Drug main - Drug - Protein main
double getdistmain(vector <double> vecDH, vector <double> vecPM, int ligid, int dm_id, int dho_id, int pm_oid, int pm_aid, double x, double y, double z)
   //d, P, d_l_id, dmainid dsideid, pmainoid, x, y, z 
   {
    double dcal,dist,xcomp,ycomp,zcomp,xc,yc,zc,ccount=0;
    dcal = 0.0;
    dist = 0.0;
    xcomp = vecDH[0] - vecPM[0];
    ycomp = vecDH[1] - vecPM[1];
    zcomp = vecDH[2] - vecPM[2];
    if (abs(xcomp) <= x/2)             xc = abs(xcomp);
    if (abs(ycomp) <= y/2)             yc = abs(ycomp);
    if (abs(zcomp) <= z/2)             zc = abs(zcomp);
    if (abs(xcomp) > (x/2))             xc = x-xcomp;
    if (abs(ycomp) > (y/2))             yc = y-ycomp;
    if (abs(zcomp) > (z/2))             zc = z-zcomp;
    dcal = dcal + (xc * xc) + (yc * yc) + (zc * zc);
    dist = sqrt(dcal);
    if ((dist < 3.00) && (dist !=0.0))
       {
       ccount = ccount + 1;
       hdist.push_back(ligid);
       hdist.push_back(dm_id);  //dmain id
       hdist.push_back(dho_id);  //dhox id
       hdist.push_back(pm_oid);
       hdist.push_back(pm_aid);
       hdist.push_back(dist);  // dist
       hbdata.push_back(hdist);
       hdist.clear();
       }
   return ccount;
   }
*/
//vec cross
vector <double> veccross(vector <double> vecA, vector <double> vecB)
   {
    vector <double> vcross;
    vcross.resize(3);
    double det1=0.0,det2=0.0,det3=0.0;
    det1 = (vecA[1]*vecB[2])-(vecA[2]*vecB[1]);
    det2 = -((vecA[0]*vecB[2])-(vecA[2]*vecB[0]));
    det3 = (vecA[0]*vecB[1])-(vecA[1]*vecB[0]);
    vcross[0] = det1;
    vcross[1] = det2;
    vcross[2] = det3;
    return vcross;
    }

//vec mag
double vecmag(vector <double> vecA)
    {
     double vmag,mag;
     vmag=0.0;
     for (int j=0; j<3; j++)         vmag = vmag + (vecA[j] * vecA[j]);
     mag = sqrt(vmag);
     return mag;
    }

vector <double> vec_norm(vector <double> vecAB)
       {
       double mag,magAB=0.0;
       mag=0.0;
       vector <double> vnorm;
       vnorm.resize(3);
       for (int k=0; k<3; k++)       mag = mag + (vecAB[k] * vecAB[k]);
       magAB=sqrt(mag);
       for (int i=0; i<vnorm.size(); i++)
           {
           vnorm[i] = vecAB[i]/magAB;
           }
       return vnorm;
       }

//vec dot
double vecdot(vector <double> vecA, vector <double> vecB)
    {
    double dot=0.0;
    for (int i=0; i<3; i++)    dot = dot + (vecA[i] * vecB[i]);
    return dot;
    } 

//vec ang
double vecang(double v_dot, double v_mag1, double v_mag2)
    {
    double ang=0.0, magn=0.0;
    magn = v_mag1 * v_mag2;
    ang = acos(v_dot/magn);
    return ang;
    }

double radtodeg(double radang)
     {
     double deg;
     deg = radang * 180 / pi;
     return deg;
     }
//Program to calculate the hydrogen bond information in DPW system
int main(int argc, char** argv)
    {
    int ligandat,proteinmainid;
    int m_ind,n_ind,n,mindex,nindex,p_orindex,p_acindex,lindex,cindex,orind,acind,pmainactual,hbonds,countwoang,cwa,flag,hbondslw,countwoanglw;
    vector <int> mainind,sideind;
    vector <double> ddhnorm,dhpnorm;
    clock_t time_req;
    float time=0.0;
    string info;
    info=argv[1];
    getinfo(info);    //info                                                   

    cout << "data" << "\t" << datain.size() << endl;
//read input info file
    nframes=atoi(datain[0].c_str());        //nf
    framestart=atoi(datain[1].c_str());     //fs
    structure=datain[2];                    //psf
    traj=datain[3];                         //dcd
    proteinid=datain[4];
    xl=atof(datain[5].c_str());
    yl=atof(datain[6].c_str());
    zl=atof(datain[7].c_str());
    proteinchain=datain[8];                 
    atindex=datain[9];                      //donor atom index file
    nligatoms=atoi(datain[10].c_str());      //nlig atoms
    nligands=atoi(datain[11].c_str());       //lnum 
    rmsd=datain[12];                         //rmsd
    mdsys=datain[13];                        //system
    inputfile=datain[14];
    inputdata=datain[15];
    //
    cout << "mains" << "\t" << mainoratomnum.size() << "\t" << sideoratomnum.size() << endl;
    cout << "nligands" << "\t" << nligands << endl;
    cout << "nf" << "\t" << nframes << endl;
    drug_index(atindex);                                            //              2
    cout << xl << "\t" << yl << "\t" << zl << endl;
    cout << "index" << "\t" << drugdonorindex.size() << "\t" << drugdonorhydindex.size() << endl;    
    //
    ofstream outfileligandsmain[nligands];
    // Total
    ofstream total;
    ostringstream tl;
    tl << "Totalmainandsidelwise_" << selection << "_" << mdsys << "_Ligands" << rmsd << "_PBC.dat";
    string totald=tl.str();
    convert(totald,total);
    // Total
    ofstream ntot;
    ostringstream totalhbonds;
    totalhbonds << "Totalhbondsframewise.dat";
    string nt=totalhbonds.str();
    convert(nt,ntot);
     // Total
    ofstream ntotlw;
    ostringstream totalhbondslw;
    totalhbondslw << "Totalhbondslwise.dat";
    string ntlw=totalhbondslw.str();
    convert(ntlw,ntotlw);
    // Total
    ofstream ntotl;
    ostringstream totalhbondsl;
    totalhbondsl << "Totall.dat";
    string ntl=totalhbondsl.str();
    convert(ntl,ntotl);
    //
    ofstream totov;
    ostringstream ov;
    ov << "Total.dat";
    string o=ov.str();
    convert(o,totov);
    //
    cout << inputfile << "\t" << inputdata << endl;
    proteinsidemain(inputfile);
    proteinmain(inputdata);
    cout << mainoratomnum.size() << endl;
    cout << proteinmainind.size() << "\t" << proteinsideind.size() << endl;
    int maintot,sidetot;
    //Total Hbonds
    //Total HBonds fo main atoms
    ofstream mstotal;
    ostringstream mst;
    mst << "Ligands_TotalHbonds_for_mainatoms_" << mdsys << "allframesrmsd" << rmsd << ".dat";
    string mstt=mst.str();
    convert(mstt,mstotal);  //
    //Total HBonds for side atoms
/*    ofstream sidetotal;
    ostringstream st;
    st << "Ligands_TotalHbonds_for_sideatoms_" << mdsys << "allframesrmsd" << rmsd << ".dat";
    string stt=st.str();
    convert(stt,sidetotal); */
/*    ofstream indexdatmain;
    ostringstream osm;
    osm << "InddatmaintoLigside" << selection << mdsys << "_rmsd" << rmsd << "_PBC.dat";
    string outr=osm.str();
    convert(outr,indexdatmain);
    ofstream indexdatside;
    ostringstream osss;
    osss << "InddatsidetoLigside" << selection << mdsys << "_rmsd" << rmsd << "_PBC.dat";
    string outs=osss.str();
    convert(outs,indexdatside);*/
    ofstream Listm;
    ostringstream List_m;
    List_m << "listmain" << selection << mdsys << "_rmsd" << rmsd << "_PBC.dat";
    string ls_m=List_m.str();
    convert(ls_m,Listm);
    ofstream Lists;
    ostringstream List_s;
    List_s << "listside" << selection << mdsys << "_rmsd" << rmsd << "_PBC.dat";
    string ls_s=List_s.str();
    convert(ls_s,Lists);
/*    int l=outr.length();
    char file[l];
    strcpy(file,outr.c_str());
    fout.open(file);*/

    ligandat = nligatoms * nligands;    //buffer atoms value
    cout << "totlig" << "\t" << ligandat << endl;
    Totalmain.resize(nligands);
    Totalside.resize(nligands);
    totalh.resize(nligands); 
    countwoangl.resize(nligands);
    for (int l=0; l<nligands; l++)    
        {
        totalh[l]=0;
        countwoangl[l]=0;
        }
    // framewise data    
    for (int i=framestart; i<=nframes; i++)          //nframes
        {
        hbonds=0; //ligand wise
	countwoang=0;
        cout << "frame" << i << endl;
        time_req = clock();
        catdcd(structure,traj,i);             //frames                              3
        ostringstream oss;
        oss << "frame.pdb";
        string filename = oss.str();

        //get protein coord
        protein_getcoord(filename,proteinchain);           //protein                4
        cout << "protein" << "\t" << procoord.size() << endl;

	//get drug coord
        drug_getcoord(filename);   //drug                                           5
        cout << "drug" << "\t" << drugcoord.size() << endl;
        maintot=0;
	sidetot=0;
        for (int j=0; j<nligands; j++)
          {
          hbdata.clear();
          hbondslw = 0;
 	  countwoanglw = 0;
          if ((i==framestart))
             {
             ostringstream ossligands1;
             ossligands1 << "Lig_" << j << "_Mainandsideallatoms_toprotein_Ligand" << selection << "_" << mdsys << "frames_rmsd" << rmsd <<"_PBC.dat";
             string outligands=ossligands1.str();
             convert(outligands,outfileligandsmain[j]);
             } //frame1 - outfiles
          for (int k=0; k<drugdonorindex.size(); k++)
              {
              m_ind = drugdonorhydindex[k] + (nligatoms * j);  //dh/o
              n_ind = drugdonorindex[k] + (nligatoms * j); // dm
              //ligand side to protein main
              for (int p=0; p<procoord.size(); p++)
                 {
                 proteinmainid = indx[p] - 1 - ligandat;                  
                 cwa = getdist_psidem_druglig(drugcoord[m_ind-1],procoord[proteinmainid],j,n_ind,m_ind,indx[p]-1,proteinmainid,xl,yl,zl);  //protein and drug donor hydrogen , oxygen
                 countwoang = countwoang + cwa;
                 }
              }
              for (int l=0; l<hbdata.size(); l++)
                  {
                  lindex=hbdata[l][0]; //drug ligand index
                  mindex=hbdata[l][1]; //drug main index
                  nindex=hbdata[l][2]; //drug ho index
                  p_orindex=hbdata[l][3]; //pmain original id  - pdb file - coordinates
 	          p_acindex=hbdata[l][4]; //pmain actual id - vector - coordinates
		  orind = p_orindex + 1; 
		  acind = p_acindex + 1;
//                  if (count(proteinmainind.begin(), proteinmainind.end(), p_acindex))  flag = 0;		    
//                  if (count(proteinsideind.begin(), proteinsideind.end(), p_acindex))  flag = 1;
                  vcddh = veccross(drugcoord[mindex-1],drugcoord[nindex-1]);
                  vcdhp = veccross(procoord[p_acindex],drugcoord[nindex-1]);
                  ddhnorm = vec_norm(vcddh);
                  dhpnorm = vec_norm(vcdhp);
                  dotp = vecdot(ddhnorm,dhpnorm);
                  magddh = vecmag(ddhnorm);
                  magdhp = vecmag(dhpnorm);
                  radangle = vecang(dotp,magddh,magdhp);
                  angle = radtodeg(radangle);
                  countwoangl[j] = countwoangl[j] + 1;
                  if (angle >= 90)
                     {
                     if (count(proteinmainind.begin(), proteinmainind.end(), acind))  flag = 0;
		     else flag = 1;
                     totalh[j] = totalh[j] + 1;
                     //if (flag ==0) cout << flag << "\t" << atoms[p_acindex] << "\t" << acind << endl;
                     if (flag == 0)
                        {
                        hbondslw = hbondslw + 1;
                        hbonds = hbonds + 1;    // total frame - wise
                        maintot = maintot + 1;  // for total main - frame wise
		        Listm << p_orindex << "\t" << lindex << endl;   // index of ligands & List of main atoms
                        Totalmain[lindex] = Totalmain[lindex] + 1;  // Total main - ligands wise
                        outfileligandsmain[lindex] << "Main" << "\t" << i << "\t" << lindex << "\t" << mindex << "\t" << dname[mindex-1] << "\t" << nindex << "\t" << dname[nindex-1] << "\t" << "Res " << res[p_acindex] << "\t" << proresid[p_acindex] << "\t" << p_orindex << "\t" << hbdata[l][5] << "\t" << angle << "\t" << acind << endl; //main angle
                        } // flag 0
                     if (flag == 1)
                        {
			hbonds = hbonds + 1;   //total frame - wise
 			sidetot = sidetot + 1; //total side - frame wise
			Lists << p_orindex << "\t" << lindex << endl;    //index of ligands & List of side atoms
                        Totalside[lindex] = Totalside[lindex] + 1; // Total side - ligands wise
                        outfileligandsmain[lindex] << "Side" << "\t" << i << "\t" << lindex << "\t" << mindex << "\t" << dname[mindex-1] << "\t" << nindex << "\t" << dname[nindex-1] << "\t" << atoms[p_acindex] << "\t" << "Res " << "\t" << proresid[p_acindex] << "\t" << res[p_acindex] << "\t" << p_orindex << "\t" << hbdata[l][5] << "\t" << angle << "\t" << acind << endl; //side angle 
                        } // flag 1
                     } // main and side
                  } // hbdata
        countwoanglw = hbdata.size();
	ntotlw << i << "\t" << countwoanglw << "\t" << hbondslw << endl;
        hbdata.clear();
        } // nligands
        mstotal << i  << "\t" << maintot << "\t" << sidetot << endl;
        procoord.clear();
        drugcoord.clear();
        hbdata.clear();
        dname.clear();
        proresid.clear();
	res.clear();    
	ntot << i << "\t" << countwoang << "\t" << hbonds << endl;
        time_req = clock() - time_req;
        time = time + ((float)time_req/CLOCKS_PER_SEC);
        }  //nframes
    // with respect to ligands
    for (int i=0; i<nligands; i++)          
        {
        ntotl << i << "\t" << countwoangl[i] << endl;    //count with only distance parameter
        total << i << "\t" << Totalmain[i] << "\t" << Totalside[i] << endl;    // count with both distance and angle parameter
	totov << totalh[i] << endl;
        } 
/*        if (selection == "main")
           {
           ligout << "Ligandmain " << i << " to proteiin main" << "\t" << Totalmain[i] << endl;
           }
        if (selection == "side")
           {
           ligout << "Ligandmain " << i << " to protein side" << "\t" << Totalside[i] << endl;
           ligoutside << "Ligandside " << i << " to protein side" << "\t" << Totalsidetoside[i] << endl;
           }*/
    cout << "Time required for " << nframes << " is " << time << " seconds" << endl;        
    return 0;
    }
