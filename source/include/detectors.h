//detector background integrated between two energies in events/t/year
double intBgRate(detector det, double Er_min, double Er_max);         

//differential background rate /keV/t/year
double diffBgRate(detector det, double Er);                            

//for declaring a new detector
int newDetector(detector *det, char *name, double exp, int ndet);
