#include "hhpred_hmm.h"

//------- amino acid utility -------//
// const char aa[]="ARNDCQEGHILKMFPSTWYVX-";
//Amino acids Sorted by alphabet     -> internal numbers a
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  X
extern const int s2a[]={ 0, 4, 3, 6,13, 7, 8, 9,11,10,12, 2,14, 5, 1,15,16,19,17,18,20};
//Internal numbers a for amino acids -> amino acids Sorted by alphabet:
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
extern const int a2s[]={ 0,14,11, 2, 1,13, 3, 5, 6, 7, 9, 8,10, 4,12,15,16,18,19,17,20};



//-------------- header definition ----------//
float Gonnet[NAA][NAA]={
//  A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V
 {10227, 3430, 2875, 3869, 1625, 2393, 4590, 6500, 2352, 3225, 5819, 4172, 1435, 1579, 3728, 4610, 6264,  418, 1824, 5709}, // A
 { 3430, 7780, 2209, 2589,  584, 2369, 3368, 3080, 2173, 1493, 3093, 5701,  763,  859, 1893, 2287, 3487,  444, 1338, 2356}, // R
 { 2875, 2209, 3868, 3601,  501, 1541, 2956, 3325, 1951, 1065, 2012, 2879,  532,  688, 1480, 2304, 3204,  219, 1148, 1759}, // N
 { 3869, 2589, 3601, 8618,  488, 2172, 6021, 4176, 2184, 1139, 2151, 3616,  595,  670, 2086, 2828, 3843,  204, 1119, 2015}, // D
 { 1625,  584,  501,  488, 5034,  355,  566,  900,  516,  741, 1336,  591,  337,  549,  419,  901, 1197,  187,  664, 1373}, // C
 { 2393, 2369, 1541, 2172,  355, 1987, 2891, 1959, 1587, 1066, 2260, 2751,  570,  628, 1415, 1595, 2323,  219,  871, 1682}, // Q
 { 4590, 3368, 2956, 6021,  566, 2891, 8201, 3758, 2418, 1624, 3140, 4704,  830,  852, 2418, 2923, 4159,  278, 1268, 2809}, // E
 { 6500, 3080, 3325, 4176,  900, 1959, 3758,26066, 2016, 1354, 2741, 3496,  741,  797, 2369, 3863, 4169,  375, 1186, 2569}, // G
 { 2352, 2173, 1951, 2184,  516, 1587, 2418, 2016, 5409, 1123, 2380, 2524,  600, 1259, 1298, 1642, 2446,  383,  876, 1691}, // H
 { 3225, 1493, 1065, 1139,  741, 1066, 1624, 1354, 1123, 6417, 9630, 1858, 1975, 2225, 1260, 1558, 3131,  417, 1697, 7504}, // I
 { 5819, 3093, 2012, 2151, 1336, 2260, 3140, 2741, 2380, 9630,25113, 3677, 4187, 5540, 2670, 2876, 5272, 1063, 3945,11005}, // L
 { 4172, 5701, 2879, 3616,  591, 2751, 4704, 3496, 2524, 1858, 3677, 7430,  949,  975, 2355, 2847, 4340,  333, 1451, 2932}, // K
 { 1435,  763,  532,  595,  337,  570,  830,  741,  600, 1975, 4187,  949, 1300, 1111,  573,  743, 1361,  218,  828, 2310}, // M
 { 1579,  859,  688,  670,  549,  628,  852,  797, 1259, 2225, 5540,  975, 1111, 6126,  661,  856, 1498, 1000, 4464, 2602}, // F
 { 3728, 1893, 1480, 2086,  419, 1415, 2418, 2369, 1298, 1260, 2670, 2355,  573,  661,11834, 2320, 3300,  179,  876, 2179}, // P
 { 4610, 2287, 2304, 2828,  901, 1595, 2923, 3863, 1642, 1558, 2876, 2847,  743,  856, 2320, 3611, 4686,  272, 1188, 2695}, // S
 { 6264, 3487, 3204, 3843, 1197, 2323, 4159, 4169, 2446, 3131, 5272, 4340, 1361, 1498, 3300, 4686, 8995,  397, 1812, 5172}, // T
 {  418,  444,  219,  204,  187,  219,  278,  375,  383,  417, 1063,  333,  218, 1000,  179,  272,  397, 4101, 1266,  499}, // W
 { 1824, 1338, 1148, 1119,  664,  871, 1268, 1186,  876, 1697, 3945, 1451,  828, 4464,  876, 1188, 1812, 1266, 9380, 2227}, // Y
 { 5709, 2356, 1759, 2015, 1373, 1682, 2809, 2569, 1691, 7504,11005, 2932, 2310, 2602, 2179, 2695, 5172,  499, 2227,11569}};// V

float Gonnet_R[NAA][NAA]={
{0.133435, 0.066867, 0.071665, 0.071669, 0.086143, 0.073328, 0.076789, 0.086161, 0.063864, 0.063859, 0.058242, 0.070022, 0.065352, 0.045193, 0.082272, 0.098916, 0.088156, 0.033515, 0.046262, 0.078575},
{0.044752, 0.151669, 0.055064, 0.047959, 0.030958, 0.072593, 0.056346, 0.040827, 0.059004, 0.029563, 0.030958, 0.095685, 0.034748, 0.024586, 0.041776, 0.049072, 0.049074, 0.035600, 0.033935, 0.032426},
{0.037511, 0.043064, 0.096418, 0.066705, 0.026559, 0.047221, 0.049453, 0.044075, 0.052976, 0.021088, 0.020138, 0.048321, 0.024228, 0.019691, 0.032662, 0.049437, 0.045091, 0.017559, 0.029116, 0.024210},
{0.050480, 0.050472, 0.089762, 0.159640, 0.025869, 0.066556, 0.100729, 0.055355, 0.059303, 0.022554, 0.021529, 0.060690, 0.027097, 0.019176, 0.046035, 0.060680, 0.054084, 0.016357, 0.028381, 0.027733},
{0.021202, 0.011385, 0.012488, 0.009040, 0.266857, 0.010878, 0.009469, 0.011930, 0.014011, 0.014673, 0.013372, 0.009919, 0.015347, 0.015713, 0.009247, 0.019333, 0.016846, 0.014994, 0.016841, 0.018897},
{0.031222, 0.046183, 0.038413, 0.040234, 0.018819, 0.060887, 0.048366, 0.025968, 0.043092, 0.021108, 0.022620, 0.046172, 0.025959, 0.017974, 0.031227, 0.034224, 0.032693, 0.017559, 0.022091, 0.023150},
{0.059887, 0.065658, 0.073684, 0.111533, 0.030004, 0.088589, 0.137200, 0.049814, 0.065657, 0.032157, 0.031428, 0.078951, 0.037799, 0.024385, 0.053362, 0.062719, 0.058531, 0.022290, 0.032160, 0.038661},
{0.084808, 0.060044, 0.082883, 0.077356, 0.047710, 0.060029, 0.062870, 0.345520, 0.054741, 0.026811, 0.027435, 0.058676, 0.033746, 0.022811, 0.052281, 0.082888, 0.058672, 0.030067, 0.030080, 0.035358},
{0.030687, 0.042362, 0.048633, 0.040456, 0.027354, 0.048630, 0.040452, 0.026723, 0.146872, 0.022237, 0.023821, 0.042362, 0.027325, 0.036034, 0.028645, 0.035232, 0.034424, 0.030709, 0.022218, 0.023274},
{0.042078, 0.029106, 0.026547, 0.021099, 0.039281, 0.032665, 0.027169, 0.017948, 0.030493, 0.127064, 0.096387, 0.031184, 0.089944, 0.063682, 0.027807, 0.033430, 0.044064, 0.033435, 0.043040, 0.103280},
{0.075922, 0.060297, 0.050153, 0.039845, 0.070823, 0.069253, 0.052531, 0.036334, 0.064625, 0.190686, 0.251356, 0.061714, 0.190682, 0.158562, 0.058923, 0.061710, 0.074195, 0.085231, 0.100056, 0.151465},
{0.054433, 0.111139, 0.071765, 0.066983, 0.031330, 0.084299, 0.078696, 0.046341, 0.068535, 0.036791, 0.036803, 0.124704, 0.043219, 0.027906, 0.051972, 0.061088, 0.061079, 0.026700, 0.036801, 0.040354},
{0.018723, 0.014874, 0.013261, 0.011022, 0.017865, 0.017466, 0.013886, 0.009822, 0.016292, 0.039107, 0.041908, 0.015928, 0.059204, 0.031798, 0.012645, 0.015942, 0.019154, 0.017479, 0.021000, 0.031793},
{0.020602, 0.016746, 0.017150, 0.012411, 0.029103, 0.019244, 0.014254, 0.010565, 0.034186, 0.044058, 0.055450, 0.016364, 0.050597, 0.175334, 0.014587, 0.018367, 0.021082, 0.080180, 0.113219, 0.035812},
{0.048640, 0.036903, 0.036892, 0.038641, 0.022212, 0.043360, 0.040452, 0.031402, 0.035245, 0.024950, 0.026724, 0.039526, 0.026095, 0.018919, 0.261161, 0.049780, 0.046442, 0.014352, 0.022218, 0.029990},
{0.060148, 0.044584, 0.057432, 0.052386, 0.047763, 0.048875, 0.048901, 0.051206, 0.044586, 0.030850, 0.028786, 0.047784, 0.033837, 0.024500, 0.051199, 0.077481, 0.065948, 0.021809, 0.030131, 0.037092},
{0.081729, 0.067978, 0.079866, 0.071188, 0.063454, 0.071183, 0.069579, 0.055262, 0.066417, 0.061998, 0.052767, 0.072842, 0.061982, 0.042875, 0.072827, 0.100547, 0.126590, 0.031831, 0.045957, 0.071184},
{0.005454, 0.008656, 0.005459, 0.003779, 0.009913, 0.006711, 0.004651, 0.004971, 0.010400, 0.008257, 0.010640, 0.005589, 0.009928, 0.028621, 0.003950, 0.005836, 0.005587, 0.328817, 0.032109, 0.006868},
{0.023798, 0.026084, 0.028616, 0.020728, 0.035199, 0.026690, 0.021213, 0.015721, 0.023786, 0.033603, 0.039486, 0.024353, 0.037708, 0.127766, 0.019332, 0.025491, 0.025501, 0.101507, 0.237902, 0.030651},
{0.074487, 0.045930, 0.043847, 0.037326, 0.072784, 0.051541, 0.046994, 0.034054, 0.045916, 0.148588, 0.110149, 0.049210, 0.105201, 0.074473, 0.048088, 0.057826, 0.072788, 0.040010, 0.056483, 0.159228}};



/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
HMM::HMM(int maxres_,int v_)
{
  v=v_;
  maxres=maxres_;
  seq = new char[maxres_];         // residues of stored sequences (first at pos 1!)
  Neff_M = new float[maxres_];     // Neff_M[i] = diversity of subalignment of seqs that have residue in col i
  Neff_I = new float[maxres_];     // Neff_I[i] = diversity of subalignment of seqs that have insert in col i
  Neff_D = new float[maxres_];     // Neff_D[i] = diversity of subalignment of seqs that have delete in col i
  name = new char[maxres_];        // Full name of first sequence of original alignment (NAME field)
  l = new int[maxres_];            // l[i] = pos. of j'th match state in aligment
  f = new float*[maxres_];         // f[i][a] = prob of finding amino acid a in column i WITHOUT pseudocounts
  g = new float*[maxres_];         // f[i][a] = prob of finding amino acid a in column i WITH pseudocounts
  p = new float*[maxres_];         // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
  tr = new float*[maxres_];        // log2 of transition probabilities M2M M2I M2D I2M I2I D2M D2D
  trr = new float*[maxres_];       // tr in lin space //__130101__//by WS
  for (int i=0; i<maxres_; i++) f[i]=new(float[NAA+3]);
  for (int i=0; i<maxres_; i++) g[i]=new(float[NAA]);
  for (int i=0; i<maxres_; i++) p[i]=new(float[NAA]);  // align memory on 16b boundaries for SSE2
  for (int i=0; i<maxres_; i++) tr[i]=new(float[NTRANS]);
  for (int i=0; i<maxres_; i++) trr[i]=new(float[NTRANS]);
  name[0]='\0';
  has_pseudocounts=false;
}

/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
HMM::~HMM()
{
  delete[] seq;
  delete[] Neff_M;
  delete[] Neff_D;
  delete[] Neff_I;
  delete[] name;
  delete[] l;
  for (int i=0; i<maxres; i++) if (f[i]) delete[] f[i]; else break;
  for (int i=0; i<maxres; i++) if (g[i]) delete[] g[i]; else break;
  for (int i=0; i<maxres; i++) if (p[i]) delete[] p[i];  else break;
  for (int i=0; i<maxres; i++) if (tr[i]) delete[] tr[i]; else break;
  for (int i=0; i<maxres; i++) if (trr[i]) delete[] trr[i]; else break;
  delete[] f;
  delete[] g;
  delete[] p;
  delete[] tr;
  delete[] trr;
}

/////////////////////////////////////////////////////////////////////////////////////
// Deep-copy constructor
/////////////////////////////////////////////////////////////////////////////////////
HMM& HMM::operator=(HMM& q)
{
  v=q.v;
  L=q.L;
  for (int i=0; i<=L+1; ++i)
    {
      for (int a=0; a<NAA; ++a)
        {
          f[i][a]=q.f[i][a];
          g[i][a]=q.g[i][a];
          p[i][a]=q.p[i][a];
        }
      for (int a=0; a<NTRANS; ++a)
	{
          tr[i][a]=q.tr[i][a];
	  trr[i][a]=q.trr[i][a];
	}
      l[i]=q.l[i];
    }
	for (int a=0; a<NAA; ++a) pb[a]=q.pb[a];
	for (int a=0; a<NAA; ++a) pav[a]=q.pav[a];
  for (int i=0; i<=L+1; ++i) Neff_M[i]=q.Neff_M[i];
  for (int i=0; i<=L+1; ++i) Neff_I[i]=q.Neff_I[i];
  for (int i=0; i<=L+1; ++i) Neff_D[i]=q.Neff_D[i];
  Neff_HMM=q.Neff_HMM;
  strcpy(name,q.name);
  strcpy(seq,q.seq);
  has_pseudocounts=q.has_pseudocounts;
  return (HMM&) (*this);
}

/////////////////////////////////////////////////////////////////////////////////////
//// Read an HMM from an HHsearch .hhm file; return 0 at end of file
/////////////////////////////////////////////////////////////////////////////////////
int HMM::Read(FILE* dbf)
{
  char line[LINELEN]="";    // input line
  char str3[8]="",str4[8]=""; // first 3 and 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  int warn=0;

	//-- INIT --//
  L=0;
  Neff_HMM=0;
  name[0]='\0';
  has_pseudocounts=false;
  //If at the end of while-loop L is still 0 then we have reached end of db file

  while (fgetline(line,LINELEN-1,dbf) && !(line[0]=='/' && line[1]=='/'))
    {

      if (strscn(line)==NULL) continue;    // skip lines that contain only white space
      substr(str3,line,0,2);               // copy the first three characters into str3
      substr(str4,line,0,3);               // copy the first four characters into str4

      /////////////////////////////////////////////////////////////////////////////////////
      // Read header for HMM
//      if (!strncmp("HH",line,2)) continue;

      if (!strcmp("NAME",str4))
        {
          ptr=strscn(line+4);              //advance to first non-white-space character
          if (ptr)
            {
              strmcpy(name,ptr,MAXRES-1); //get name
            }
          else
            {
              strcpy(name,"undefined");
            }
          if (v>=4) cout<<"Reading in HMM "<<name<<":\n";
        }
      else if (!strcmp("SEQ",str3))
      	{
      		ptr=strscn(line+4);              //advance to first non-white-space character
          if (ptr)
            {
              strmcpy(seq,ptr,MAXRES-1); //get sequence
            }
          else
            {
              strcpy(seq,"undefined");
            }
          if (v>=4) cout<<"Reading in SEQ "<<seq<<":\n";
      	}
      else if (!strcmp("LENG",str4))
        {
          ptr=line+4;
          L=strint(ptr);                   //read next integer (number of match states)
        }
//      else if (!strcmp("NEFF",str4)) sscanf(line+6,"%f",&Neff_HMM);
//      else if (!strcmp("DESC",str4)) continue;
//      else if (!strcmp("COM",str3))  continue;
//      else if (!strcmp("DATE",str4)) continue;


      /////////////////////////////////////////////////////////////////////////////////////
      // Read average amino acid frequencies for HMM
      else if (!strcmp("NULL",str4))
        {
          ptr=line+4;
          for (a=0; a<NAA && ptr; ++a)
            //s2[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
            pb[s2a[a]] = (float) fpow2(float(-strinta(ptr))/HMMSCALE);
          if (!ptr) return Warning(dbf,line,name);
          if (v>=4)
            {
              printf("\nNULL  ");
              for (a=0; a<NAA; ++a) printf("%5.1f ",100.*pb[s2a[a]]);
              printf("\n");
            }
        }

      /////////////////////////////////////////////////////////////////////////////////////
      // Read transition probabilities from start state
      else if (!strcmp("HMM",str3))
        {
          fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
          fgetline(line,LINELEN-1,dbf); // Skip line with transition labels
          ptr=line;

          for (a=0; a<=D2D && ptr; ++a)
            tr[0][a] = float(-strinta(ptr))/HMMSCALE; //store transition probabilites as log2 values
	  // strinta returns next integer in string and puts ptr to first char
	  // after the integer. Returns -99999 if '*' is found.
	  // ptr is set to 0 if no integer is found after ptr.
          Neff_M[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
          Neff_I[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
          Neff_D[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
          if (!ptr) return Warning(dbf,line,name);

          /////////////////////////////////////////////////////////////////////////////////////
          // Read columns of HMM
          int next_i=0;  // index of next column
          while (fgetline(line,LINELEN-2,dbf) &&  !(line[0]=='/' && line[1]=='/') && line[0]!='#')
            {
              if (strscn(line)==NULL) continue; // skip lines that contain only white space

              // Read in AA probabilities
              ptr=line+1;
              int prev_i = next_i;
              next_i = strint(ptr); ++i;
              if (v && next_i!=prev_i+1)
                if (++warn<=5)
                  {
//                    cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<prev_i<<" is followed by state "<<next_i<<"\n";
//                    if (warn==5) cerr<<endl<<"WARNING: further warnings while reading HMMs will be suppressed.\n";
                  }
              if (i>L)
                {
                  cerr<<endl<<"WARNING: in HMM "<<name<<" there are more columns than the stated length "<<L<<". Skipping HMM\n";
                  return 2;
                }
              if (i>MAXRES-2)
                {
                  fgetline(line,LINELEN-1,dbf); // Skip line
                  continue;
                }

              for (a=0; a<NAA && ptr; ++a)
                f[i][s2a[a]] = fpow2(float(-strinta(ptr))/HMMSCALE);       // speed-up ~5 s for 10000 SCOP domains

              //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
              l[i]=strint(ptr);
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("%s",line);
                  printf("%6i ",i);
                  for (a=0; a<NAA; ++a) printf("%5.1f ",100*f[i][s2a[a]]);
                  printf("%5i",l[i]);
                  printf("\n");
                }

              // Read transition probabilities
              fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
              if (line[0]!=' ' && line[0]!='\t') return Warning(dbf,line,name);
              ptr=line;
              for (a=0; a<=D2D && ptr; ++a)
                tr[i][a] = float(-strinta(ptr))/HMMSCALE; //store transition prob's as log2-values
              Neff_M[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
              if (Neff_M[i] == 0) { Neff_M[i] = 1; }
              Neff_I[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
              Neff_D[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
		Neff_HMM+=Neff_M[i];
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("       ");
                  for (a=0; a<=D2D; ++a) printf("%5.1f ",100*fpow2(tr[i][a]));
                  printf("%5.1f %5.1f %5.1f \n",Neff_M[i],Neff_I[i],Neff_D[i]);
                }
            }
          if (line[0]=='/' && line[1]=='/') break;
        }
      else if (v) cerr<<endl<<"WARNING: Ignoring line\n\'"<<line<<"\'\nin HMM "<<name<<"\n";

    } //while(getline)


  if (L==0) return 0; //End of db file -> stop reading in
	Neff_HMM/=L;

  if (v && (int)strlen(seq)!=L) {cerr<<endl<<"ERROR: in HMM "<<name<<" length not equal "<<strlen(seq)<<" != "<<L<<"\n";return -1;}
  if (v && i!=L) {cerr<<endl<<"WARNING: in HMM "<<name<<" there are only "<<i<<" columns while the stated length is "<<L<<"\n";return -1;}
  if (v && i>MAXRES-2) {i=MAXRES-2; cerr<<endl<<"WARNING: maximum number "<<MAXRES-2<<" of residues exceeded while reading HMM "<<name<<"\n";return -1;}
  if (v && !i) {cerr<<endl<<"WARNING: HMM "<<name<<" contains no match states. Check the alignment that gave rise to this HMM.\n";return -1;}
  if (v>=2) {cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";}
  L = i;
  
// Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<NAA; ++a) f[0][a]=f[L+1][a]=pb[a];
  Neff_M[L+1]=1.0f;
  Neff_I[L+1]=Neff_D[L+1]=0.0f;

  return 1; //return status: ok
}

//================ WS_Additional_Read =============//__130131__//
int HMM::Read_HMM(FILE* dbf,int L)
{
  char line[LINELEN]="";    // input line
  char str3[8]="",str4[8]=""; // first 3 and 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  int warn=0;

	//-- INIT --//
  Neff_HMM=0;
//  name[0]='\0';
  has_pseudocounts=false;

	//---- process null -----//
	fgetline(line,LINELEN-1,dbf);
	fgetline(line,LINELEN-1,dbf); //-> NULL
	{
    ptr=line+4;
    for (a=0; a<NAA && ptr; ++a)
      //s2[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
      pb[s2a[a]] = (float) fpow2(float(-strinta(ptr))/HMMSCALE);
    if (!ptr) return Warning(dbf,line,name);
    if (v>=4)
      {
        printf("\nNULL  ");
        for (a=0; a<NAA; ++a) printf("%5.1f ",100.*pb[s2a[a]]);
        printf("\n");
      }
	}
	fgetline(line,LINELEN-1,dbf);  //-> HMM

	/////////////////////////////////////////////////////////////////////////////////////
	// Read transition probabilities from start state
	  fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
	  fgetline(line,LINELEN-1,dbf); // Skip line with transition labels
	  ptr=line;
	  for (a=0; a<=D2D && ptr; ++a)
	    tr[0][a] = float(-strinta(ptr))/HMMSCALE; //store transition probabilites as log2 values
// strinta returns next integer in string and puts ptr to first char
// after the integer. Returns -99999 if '*' is found.
// ptr is set to 0 if no integer is found after ptr.
	  Neff_M[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
	  Neff_I[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
	  Neff_D[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
	  if (!ptr) return Warning(dbf,line,name);
	
	  /////////////////////////////////////////////////////////////////////////////////////
	  // Read columns of HMM
	  int next_i=0;  // index of next column
	int count=0;
	  while (fgetline(line,LINELEN-2,dbf) &&  !(line[0]=='/' && line[1]=='/') && line[0]!='#')
	    {
		count++;
		if(count>3*L)break;
	      if (strscn(line)==NULL) continue; // skip lines that contain only white space
	
	      // Read in AA probabilities
	      ptr=line+1;
	      int prev_i = next_i;
	      next_i = strint(ptr); ++i;
	      if (v && next_i!=prev_i+1)
	        if (++warn<=5)
	          {
//	            cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<prev_i<<" is followed by state "<<next_i<<"\n";
//	            if (warn==5) cerr<<endl<<"WARNING: further warnings while reading HMMs will be suppressed.\n";
	          }
	      if (i>L)
	        {
	          cerr<<endl<<"WARNING: in HMM "<<name<<" there are more columns than the stated length "<<L<<". Skipping HMM\n";
	          return 2;
	        }
	      if (i>MAXRES-2)
	        {
	          fgetline(line,LINELEN-1,dbf); // Skip line
	          continue;
	        }
	
	      for (a=0; a<NAA && ptr; ++a)
	        f[i][s2a[a]] = fpow2(float(-strinta(ptr))/HMMSCALE);       // speed-up ~5 s for 10000 SCOP domains
	
	      //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
	      l[i]=strint(ptr);
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4)
	        {
	          printf("%s",line);
	          printf("%6i ",i);
	          for (a=0; a<NAA; ++a) printf("%5.1f ",100*f[i][s2a[a]]);
	          printf("%5i",l[i]);
	          printf("\n");
	        }
	
	      // Read transition probabilities
	      fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
	      if (line[0]!=' ' && line[0]!='\t') return Warning(dbf,line,name);
	      ptr=line;
	      for (a=0; a<=D2D && ptr; ++a)
	        tr[i][a] = float(-strinta(ptr))/HMMSCALE; //store transition prob's as log2-values
	      Neff_M[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
	      if (Neff_M[i] == 0) { Neff_M[i] = 1; }
	      Neff_I[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
	      Neff_D[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
	      Neff_HMM+=Neff_M[i];
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4)
	        {
	          printf("       ");
	          for (a=0; a<=D2D; ++a) printf("%5.1f ",100*fpow2(tr[i][a]));
	          printf("%5.1f %5.1f %5.1f \n",Neff_M[i],Neff_I[i],Neff_D[i]);
	        }
	    }
	 
	//----- final process ---//
	Neff_HMM/=L;
  if (v && (int)strlen(seq)!=L) {cerr<<endl<<"ERROR: in HMM "<<name<<" length not equal "<<strlen(seq)<<" != "<<L<<"\n";return -1;}
  if (v && i!=L) {cerr<<endl<<"WARNING: in HMM "<<name<<" there are only "<<i<<" columns while the stated length is "<<L<<"\n";return -1;}
  if (v && i>MAXRES-2) {i=MAXRES-2; cerr<<endl<<"WARNING: maximum number "<<MAXRES-2<<" of residues exceeded while reading HMM "<<name<<"\n";return -1;}
  if (v && !i) {cerr<<endl<<"WARNING: HMM "<<name<<" contains no match states. Check the alignment that gave rise to this HMM.\n";return -1;}
  if (v>=2) {cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";}
  
// Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<NAA; ++a) f[0][a]=f[L+1][a]=pb[a];
  Neff_M[L+1]=1.0f;
  Neff_I[L+1]=Neff_D[L+1]=0.0f;
  return 1; //return status: ok
}

int HMM::Read_TPL(FILE* dbf)
{
  char line[LINELEN]="";    // input line
  char str3[255]="",str4[255]=""; // first 3 and 4 letters of input line
  char strx[255]="";
  char* ptr;                // pointer for string manipulation

  while (fgetline(line,LINELEN-1,dbf) )
    {

      if (strscn(line)==NULL) continue;    // skip lines that contain only white space
      substr(str3,line,0,12);               // copy the first three characters into str3
      substr(str4,line,0,5);               // copy the first four characters into str4
      substr(strx,line,0,29);

      /////////////////////////////////////////////////////////////////////////////////////
      // Read header for HMM

      if (!strcmp("Template Name",str3))
        {
          ptr=strscn(line+17);              //advance to first non-white-space character
          if (ptr)
            {
              strmcpy(name,ptr,MAXRES-1); //get name
            }
          else
            {
              strcpy(name,"undefined");
            }
          if (v>=4) cout<<"Reading in HMM "<<name<<":\n";
        }
      else if (!strcmp("SEQRES sequen",str3))
      	{
      		ptr=strscn(line+18);              //advance to first non-white-space character
          if (ptr)
            {
              strmcpy(seq,ptr,MAXRES-1); //get sequence
            }
          else
            {
              strcpy(seq,"undefined");
            }
          if (v>=4) cout<<"Reading in SEQ "<<seq<<":\n";
      	}
      else if (!strcmp("Length",str4))
        {
          ptr=line+11;
          L=strint(ptr);                   //read next integer (number of match states)
        }
//      else if (!strcmp("NEFF",str4)) sscanf(line+6,"%f",&Neff_HMM);
//      else if (!strcmp("DESC",str4)) continue;
//      else if (!strcmp("COM",str3))  continue;
//      else if (!strcmp("DATE",str4)) continue;

			//----- check HMM -----//
			else if(!strcmp("//////////// Original HHM file",strx))
			{
				Read_HMM(dbf,L);
				return 1;
			}
	}
  return 0; //return status: ok
}

int HMM::Read_TGT(FILE* dbf)
{
  char line[LINELEN]="";    // input line
  char str3[255]="",str4[255]=""; // first 3 and 4 letters of input line
  char strx[255]="";
  char* ptr;                // pointer for string manipulation

  while (fgetline(line,LINELEN-1,dbf) )
    {
      if (strscn(line)==NULL) continue;    // skip lines that contain only white space
      substr(str3,line,0,12);               // copy the first three characters into str3
      substr(str4,line,0,5);               // copy the first four characters into str4
      substr(strx,line,0,29);


      /////////////////////////////////////////////////////////////////////////////////////
      // Read header for HMM

      if (!strcmp("Sequence Name",str3))
        {
          ptr=strscn(line+17);              //advance to first non-white-space character
          if (ptr)
            {
              strmcpy(name,ptr,MAXRES-1); //get name
            }
          else
            {
              strcpy(name,"undefined");
            }
          if (v>=4) cout<<"Reading in HMM "<<name<<":\n";
        }
      else if (!strcmp("Sequen",str4))
      	{
      		ptr=strscn(line+11);              //advance to first non-white-space character
          if (ptr)
            {
              strmcpy(seq,ptr,MAXRES-1); //get sequence
            }
          else
            {
              strcpy(seq,"undefined");
            }
          if (v>=4) cout<<"Reading in SEQ "<<seq<<":\n";
      	}
      else if (!strcmp("Length",str4))
        {
          ptr=line+11;
          L=strint(ptr);                   //read next integer (number of match states)
        }
//      else if (!strcmp("NEFF",str4)) sscanf(line+6,"%f",&Neff_HMM);
//      else if (!strcmp("DESC",str4)) continue;
//      else if (!strcmp("COM",str3))  continue;
//      else if (!strcmp("DATE",str4)) continue;

			//----- check HMM -----//
			else if(!strcmp("//////////// Original HHM file",strx))
			{
				Read_HMM(dbf,L);
				return 1;
			}
	}
  return 0; //return status: ok
}


/////////////////////////////////////////////////////////////////////////////////////
// Add transition pseudocounts to HMM (and calculate lin-space transition probs)
/////////////////////////////////////////////////////////////////////////////////////
void HMM::AddTransitionPseudocounts(float gapd, float gape, float gapf, float gapg, float gaph, float gapi, float gapb)
{
  // init check
  if (has_pseudocounts)return;
  

  int i;               //position in alignment
  float sum;
  float pM2M, pM2I, pM2D, pI2I, pI2M, pD2D, pD2M;
  float p0,p1,p2;

  // Calculate pseudocount transition probabilities
  pM2D=pM2I=gapd*0.0286;     //a-priori probability for inserts and deletions
  pM2M=1-pM2D-pM2I;
  // gape=0 -> pI2I=0   gape=1 -> pI2I=0.75    gape=inf -> pI2I=1.
  pI2I=1.0*gape/(gape-1+1.0/0.75);
  pI2M=1-pI2I;
  // gape=0 -> pD2D=0   gape=1 -> pD2D=0.75    gape=inf -> pD2D=1.
  pD2D=1.0*gape/(gape-1+1.0/0.75);
  pD2M=1-pD2D;

  for (i=0; i<=L; ++i) //for all columns in HMM
    {
      // Transitions from M state
      p0 = (Neff_M[i]-1)*fpow2(tr[i][M2M]) + gapb*pM2M;
      p1 = (Neff_M[i]-1)*fpow2(tr[i][M2D]) + gapb*pM2D;
      p2 = (Neff_M[i]-1)*fpow2(tr[i][M2I]) + gapb*pM2I;
      if (i==0) p1=p2=0;       //from M(0) no transition to D(1) and I(0) possible
      if (i==L) p1=p2=0;       //from M(L) no transition to D(L+1) and I(L+1) possible
      sum = p0+p1+p2+FLT_MIN;
      tr[i][M2M] = fast_log2(p0/sum);
      tr[i][M2D] = fast_log2(p1/sum)*gapf;
      tr[i][M2I] = fast_log2(p2/sum)*gapg;

      // Transitions from I state
      p0 = Neff_I[i]*fpow2(tr[i][I2M]) + gapb*pI2M;
      p1 = Neff_I[i]*fpow2(tr[i][I2I]) + gapb*pI2I;
      sum = p0+p1+FLT_MIN;
      tr[i][I2M] = fast_log2(p0/sum);
      tr[i][I2I] = fast_log2(p1/sum)*gapi;

      // Transitions from D state
      p0 = Neff_D[i]*fpow2(tr[i][D2M]) + gapb*pD2M;
      p1 = Neff_D[i]*fpow2(tr[i][D2D]) + gapb*pD2D;
      if (i==L) p1=0;          //from D(L) no transition to D(L+1) possible
      sum = p0+p1+FLT_MIN;
      tr[i][D2M] = fast_log2(p0/sum);
      tr[i][D2D] = fast_log2(p1/sum)*gaph;

    }

  if (v>=4)
    {
      printf("\nPseudocount transition probabilities:\n");
      printf("pM2M=%4.1f%%, pM2I=%4.1f%%, pM2D=%4.1f%%, ",100*pM2M,100*pM2I,100*pM2D);
      printf("pI2M=%4.1f%%, pI2I=%4.1f%%, ",100*pI2M,100*pI2I);
      printf("pD2M=%4.1f%%, pD2D=%4.1f%% ",100*pD2M,100*pD2D);
      printf("tau = %4.1f%%\n\n",100.*gapb/(Neff_HMM-1+gapb));
      printf("Listing transition probabilities WITH pseudocounts:\n");
      printf("   i dssp pred sacc     M->M   M->I   M->D   I->M   I->I   D->M   D->D\n");

      for (i=1; i<=L; ++i) //for all columns in HMM
        {
		char wsc='-';
          printf("%4i  %1c    %1c    %1c    %6.3f %6.3f %6.3f ",i,wsc,wsc,wsc,fpow2(tr[i][M2M]),fpow2(tr[i][M2I]),fpow2(tr[i][M2D]));
          printf("%6.3f %6.3f ",fpow2(tr[i][I2M]),fpow2(tr[i][I2I]));
          printf("%6.3f %6.3f ",fpow2(tr[i][D2M]),fpow2(tr[i][D2D]));
	printf("\n");
//          printf("%1i %2i  %1i\n",ss_pred[i],ss_conf[i],ss_dssp[i]);
        }
      printf("\n");
//      printf("nss_dssp=%i  nss_pred=%i\n",nss_dssp,nss_pred);
    }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Generate an amino acid frequency matrix g[][] with full pseudocount admixture (tau=1)
/////////////////////////////////////////////////////////////////////////////////////
void HMM::PreparePseudocounts()
{
  for (int i=0; i<=L+1; ++i)
      for (int a=0; a<NAA; ++a)
	g[i][a] = ScalarProd20(Gonnet_R[a],f[i]);
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate amino acid background frequencies in HMM
/////////////////////////////////////////////////////////////////////////////////////
void HMM::CalculateAminoAcidBackground()
{
  int a, i;
  // initialize vector of average aa freqs with pseudocounts
  for (a=0; a<NAA; ++a) pav[a]=pb[a]*100.0f/Neff_HMM;
  // calculate averages
  for (i=1; i<=L; ++i)
    for (a=0; a<NAA; ++a)
      pav[a] += p[i][a];
  // Normalize vector of average aa frequencies pav[a]
  NormalizeTo1(pav,NAA);
  for (a=0; a<NAA; ++a) p[0][a] = p[L+1][a] = pav[a];

}

/////////////////////////////////////////////////////////////////////////////////////
// Add amino acid pseudocounts to HMM and calculate average protein aa probabilities pav[a]
// Pseudocounts: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
/////////////////////////////////////////////////////////////////////////////////////
void HMM::AddAminoAcidPseudocounts(char pcm, float pca, float pcb, float pcc)
{
  int i;               //position in HMM
  int a;               //amino acid (0..19)
  float sum;
  float tau;           //tau = pseudocount admixture

  if (has_pseudocounts) {
    pcm = 0;
  }

  // Calculate amino acid frequencies p[i][a] = (1-tau(i))*f[i][a] + tau(i)*g[i][a]
  switch (pcm)
    {
    case 0: //no pseudocounts whatsoever: tau=0
      for (i=1; i<=L; ++i)
        for (a=0; a<NAA; ++a)
          p[i][a]=f[i][a];
      break;
    case 1: //constant pseudocounts (for optimization): tau = pca
      tau = pca;
      for (i=1; i<=L; ++i)
        for (a=0; a<NAA; ++a)
          p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
      break;
    case 2: // diversity-dependent (i.e, Neff_M[i]-dependent) pseudocounts
      if (pcc==1.0f)
        for (i=1; i<=L; ++i)
          {
            tau = fmin(1.0, pca/(1. + Neff_M[i]/pcb ) );
            for (a=0; a<NAA; ++a)
              p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
          }
      else
        for (i=1; i<=L; ++i)
          {
            tau = fmin(1.0, pca/(1. + pow((Neff_M[i])/pcb,pcc)));
            for (a=0; a<NAA; ++a)
              p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
          }
      break;
    case 3: // constant-diversity pseudocounts // Is this still used? => scrap? (JS) 
      for (i=1; i<=L; ++i)
        {
          float x = Neff_M[i]/pcb;
          pca = 0.793 + 0.048*(pcb-10.0);
          tau = fmax(0.0, pca*(1-x + pcc*x*(1-x)) );
          for (a=0; a<NAA; ++a)
            p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
        }
      if (v>=2) { printf("Divergence before / after addition of amino acid pseudocounts: %5.2f / %5.2f\n",Neff_HMM, CalcNeff()); }
      break;
    } //end switch (pcm)


  //turn on pseudocount switch to indicate that HMM contains pseudocounts
  if (pcm!=0) has_pseudocounts=true;

  // DEBUGGING output
  if (v>=3)
    {
      switch (pcm)
        {
        case 0:
          cout<<"No pseudocounts added (-pcm 0)\n";
          return;
        case 1:
          cout<<"Adding constant AA pseudocount admixture of "<<pca<<" to HMM "<<name<<"\n";
          break;
        case 2:
          cout<<"Adding divergence-dependent AA pseudocounts (-pcm 2) with admixture of "
              <<fmin(1.0, pca/(1. + Neff_HMM/pcb ) )<<" to HMM "<<name<<"\n";
          break;
        } //end switch (pcm)
      if (v>=4)
        {
          cout<<"\nAmino acid frequencies WITHOUT pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
          for (i=1; i<=L; ++i)
            {
              printf("%3i:  ",i);
              sum=0;
              for (a=0; a<NAA; ++a)
                {
                  sum+=f[i][a];
                  printf("%4.1f ",100*f[i][a]);
                }
              printf("  sum=%5.3f\n",sum);
            }
          cout<<"\nAmino acid frequencies WITH pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
          for (i=1; i<=L; ++i)
            {
              printf("%3i:  ",i);
              sum=0;
              for (a=0; a<NAA; ++a)
                {
                  sum+=p[i][a];
                  printf("%4.1f ",100*p[i][a]);
                }
              printf("  sum=%5.3f\n",sum);
            }
        }
    }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Factor Null model into HMM t
// !!!!! ATTENTION!!!!!!!  after this t->p is not the same as after adding pseudocounts !!!
/////////////////////////////////////////////////////////////////////////////////////
void HMM::IncludeNullModelInHMM(HMM* q, HMM* t, int columnscore )
{

  int j;         //query and template match state indices
  int a;           //amino acid index

  // Multiply template frequencies with amino acid weights = 1/background_freq(a) (for all but SOP scores)
  switch (columnscore)
    {
    default:
    case 0: // Null model with background prob. from database
      for (j=0; j<=t->L+1; ++j)
				for (a=0; a<NAA; ++a)
	  			t->p[j][a] /= pb[a];
      break;

    case 1: // Null model with background prob. equal average from query and template
      float pnul[NAA]; // null model probabilities used in comparison (only set in template/db HMMs)
      for (a=0; a<NAA; ++a) pnul[a] = 0.5*(q->pav[a]+t->pav[a]);
      for (j=0; j<=t->L+1; ++j)
				for (a=0; a<NAA; ++a)
	  			t->p[j][a] /= pnul[a];
      break;

    case 2: // Null model with background prob. from template protein
      for (j=0; j<=t->L+1; ++j)
				for (a=0; a<NAA; ++a)
	 			 t->p[j][a] /= t->pav[a];
      break;

    case 3: // Null model with background prob. from query protein
      for (j=0; j<=t->L+1; ++j)
				for (a=0; a<NAA; ++a)
	  			t->p[j][a] /= q->pav[a];
      break;
    }

  if (v>=4)
    {
      cout<<"\nAverage amino acid frequencies\n";
      cout<<"         A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      cout<<"Q:    ";
      for (a=0; a<NAA; ++a) printf("%4.1f ",100*q->pav[a]);
      cout<<"\nT:    ";
      for (a=0; a<NAA; ++a) printf("%4.1f ",100*t->pav[a]);
      cout<<"\npb:   ";
      for (a=0; a<NAA; ++a) printf("%4.1f ",100*pb[a]);
    }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transform log to lin transition probs
/////////////////////////////////////////////////////////////////////////////////////
void HMM::Log2LinTransitionProbs(float beta)
{
  for (int i=0; i<=L; ++i)
    {
      for (int a=0; a<NTRANS; ++a)
        trr[i][a] = fpow2(beta*tr[i][a]);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Set query columns in His-tags etc to Null model distribution
/////////////////////////////////////////////////////////////////////////////////////
void HMM::NeutralizeTags()
{
  char* qseq = seq;
  char* pt;
  int a,i;

  // Neutralize His tag
  if ( (pt=strstr(qseq,"HHHHH")) )
    {
      int i0 = pt-qseq+1;
      for (i=imax(i0-8,1); i<i0; ++i)   // neutralize leading 5 columns
        for (a=0; a<NAA; ++a) p[i][a]=f[i][a]=pb[a];
      for (; (*pt)=='H'; ++i,++pt)      // neutralize His columns
        for (a=0; a<NAA; ++a) p[i][a]=f[i][a]=pb[a];
      int i1=i;
      for (; i<imin(i1+8,L+1); ++i)    // neutralize trailing 5 columns
        for (a=0; a<NAA; ++a) p[i][a]=f[i][a]=pb[a];
      if (v>=2) printf("Neutralized His-tag between positions %i and %i\n",imax(i0-8,1),i-1);
    }

  // Neutralize C-myc tag
  if ( (pt=strstr(qseq,"EQKLISEEDL")) )
    {
      if (v>=2) printf("Neutralized C-myc-tag at position %i\n",int(pt-qseq)+1);
      for (i=pt-qseq+1; i<=pt-qseq+10; ++i)
        for (a=0; a<NAA; ++a) p[i][a]=f[i][a]=pb[a];
    }
  // Neutralize FLAG tag
  if ( (pt=strstr(qseq,"DYKDDDDK")) )
    {
      if (v>=2) printf("Neutralized FLAG-tag at position %i\n",int(pt-qseq)+1);
      for (i=pt-qseq+1; i<=pt-qseq+8; ++i)
        for (a=0; a<NAA; ++a) p[i][a]=f[i][a]=pb[a];
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate effective number of sequences using profiles INCLUDING pseudocounts
/////////////////////////////////////////////////////////////////////////////////////
float HMM::CalcNeff()
{
  float Neff=0;
  for (int i=1; i<=L; ++i)
    for (int a=0; a<NAA; ++a)
      if (p[i][a]>1E-10) Neff-=p[i][a]*fast_log2(p[i][a]);
  return fpow2(Neff/L);
}

/////////////////////////////////////////////////////////////////////////////////////
// Write HMM to output file
/////////////////////////////////////////////////////////////////////////////////////
void HMM::WriteToFile(char* outfile)
{
  int i,a;
  FILE *outf=NULL;
  outf=fopen(outfile,"wb");
  if (v>=2) cout<<"Writing HMM to "<<outfile<<"\n";

  //header
  fprintf(outf,"NAME  %s\n",name);   // name 
  fprintf(outf,"LENG  %d\n",L);      // length
  fprintf(outf,"SEQ   %s\n",seq);    // sequence

  //print null model background probabilities from substitution matrix
  fprintf(outf,"NULL   ");
  for (a=0; a<NAA; ++a) fout(outf,-iround(fast_log2(pb[s2a[a]])*HMMSCALE ));
  fprintf(outf,"\n");

  // print table header line with amino acids
  fprintf(outf,"HMM    ");
  for (a=0; a<NAA; ++a) fprintf(outf,"%1c\t",i2aa(s2a[a]));
  fprintf(outf,"\n");

  // print table header line with state transitions
  fprintf(outf,"       M->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D\n");

  // print out transition probabilities from begin state (virtual match state)
  fprintf(outf,"       ");
  for (a=0; a<=D2D; ++a) fout(outf,-iround(tr[0][a]*HMMSCALE));
  fout(outf,iround(Neff_M[0]*HMMSCALE));
  fout(outf,iround(Neff_I[0]*HMMSCALE));
  fout(outf,iround(Neff_D[0]*HMMSCALE));
  fprintf(outf,"\n");

  // Start loop for printing HMM columns
  for (i=1; i<=L; ++i)
    {
      // Print amino acid and number
      fprintf(outf,"%1c %-4i ",seq[i-1],i);

      // Print emission probabilities for match state
      for (a=0; a<NAA; ++a) fout(outf,-iround(fast_log2(p[i][s2a[a]])*HMMSCALE ));
      fprintf(outf,"%-i",l[i]);
      fprintf(outf,"\n");

      // Print transition probabilities
      fprintf(outf,"       ");
      for (a=0; a<=D2D; ++a) fout(outf,-iround(tr[i][a]*HMMSCALE));
      fout(outf,iround(Neff_M[i]*HMMSCALE));
      fout(outf,iround(Neff_I[i]*HMMSCALE));
      fout(outf,iround(Neff_D[i]*HMMSCALE));
      fprintf(outf,"\n\n");
    } // end for(i)-loop for printing HMM columns

  fprintf(outf,"//\n");
  fclose(outf);
}
