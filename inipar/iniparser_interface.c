
#include <stdio.h>
#include <stdlib.h>
#include "iniparser.h"

void config_getdouble_(double* value, double *defaultValue, const char key[100], const char conf[1000]) {
	dictionary * dict = iniparser_load(conf);

	if (dict == NULL) {
		*value = *defaultValue;
	} else {
		*value = iniparser_getdouble(dict, key, *defaultValue);
		iniparser_freedict(dict);
	}
}

void config_getint_(int* value, int *defaultValue, char key[100], const char conf[1000]) {
	dictionary * dict = iniparser_load(conf);
	if (dict == NULL) {
		*value = *defaultValue;
	} else {
		*value = iniparser_getint(dict, key, *defaultValue);
		iniparser_freedict(dict);
	}
}

void config_getstr_(char value[500], char defaultValue[500], char key[1000], const char conf[1000]) {
	dictionary * dict = iniparser_load(conf);
	if (dict == NULL) {
		strncpy(value, defaultValue, 500);
	} else {
		char * str;
		str = iniparser_getstring(dict, key, "NOT VALID KEY");
   		if (str=="NOT VALID KEY") strncpy(value, defaultValue, 500);
		else strncpy(value, str, 500);
		iniparser_freedict(dict);
	}
}

void char_to_arrayint_(char str[], int *array){
	char *pt;
	int i;
	i=0;
	pt = strtok(str, ", ");

	while (pt != NULL) {
		array[i] = atoi(pt);
		pt = strtok (NULL, ", ");
		i++;
	}
}

void char_to_arraydouble_(char str[], double *array){
	char *pt;
	int i;
	i=0;
	pt = strtok(str,", ");
	while (pt != NULL) {
		array[i] = atof(pt);
		pt = strtok (NULL, ", ");
		i++;
	}
}

void config_validate_(const char conf[10000]) {
	dictionary * dict = iniparser_load(conf);

	int i, k, m, secCount, keysCount, found;
	char* secName;
	char** keys;
	const int expectedCount = 60;
	char* expected[expectedCount];

	if (dict == NULL)
		return;

	expected[0] = "mcluster:n";
	expected[1] = "mcluster:fracb";
	expected[2] = "mcluster:fracb_reference";
	expected[3] = "mcluster:initialmodel";
	expected[4] = "mcluster:w0";
	expected[5] = "mcluster:S";
	expected[6] = "mcluster:fractal";
	expected[7] = "mcluster:qvir";
	expected[8] = "mcluster:mfunc";
	expected[9] = "mcluster:single_mass";
	expected[10] = "mcluster:mlow";
	expected[11] = "mcluster:mup";
	expected[12] = "mcluster:alpha_imf";
	expected[13] = "mcluster:mlim_imf";
	expected[14] = "mcluster:alpha_l3";
	expected[15] = "mcluster:beta_l3";
	expected[16] = "mcluster:mu_l3";
	expected[17] = "mcluster:pairing";
	expected[18] = "mcluster:adis";
	expected[19] = "mcluster:eigen";
	expected[20] = "mcluster:amin";
	expected[21] = "mcluster:amax";
	expected[22] = "mcluster:tf";
	expected[23] = "mcluster:rbar";
	expected[24] = "mcluster:rh_mcl";
	expected[25] = "mcluster:conc_pop";
	expected[26] = "mcluster:potential_energy";
	expected[27] = "mcluster:epoch";
	expected[28] = "mcluster:zini";
	expected[29] = "mcluster:seedmc";
	expected[30] = "mcluster:outputf";
	expected[31] = "mcluster:check_en";
	expected[32] = "mcluster:BSE";
	expected[33] = "mcluster:pts1";
	expected[34] = "mcluster:pts2";
	expected[35] = "mcluster:pts3";
	expected[36] = "mcluster:mdflag";
	expected[37] = "mcluster:neta";
	expected[38] = "mcluster:bwind";
	expected[39] = "mcluster:hewind";
	expected[40] = "mcluster:nsflag";
	expected[41] = "mcluster:psflag";
	expected[42] = "mcluster:ecflag";
	expected[43] = "mcluster:ifflag";
	expected[44] = "mcluster:wdflag";
	expected[45] = "mcluster:mxns";
	expected[46] = "mcluster:bhflag";
	expected[47] = "mcluster:kmech";
	expected[48] = "mcluster:disp";
	expected[49] = "mcluster:bhspin";
	expected[50] = "mcluster:beta";
	expected[51] = "mcluster:xi";
	expected[52] = "mcluster:acc2";
	expected[53] = "mcluster:epsnov";
	expected[54] = "mcluster:eddfac";
	expected[55] = "mcluster:gamma";
	expected[56] = "mcluster:ceflag";
	expected[57] = "mcluster:tflag";
	expected[58] = "mcluster:alpha1";
	expected[59] = "mcluster:lambda";	

	secCount = iniparser_getnsec(dict);
	for (i = 0; i < secCount; ++i) {
		secName = iniparser_getsecname(dict, i);

		keys = iniparser_getseckeys(dict, secName);
		keysCount = iniparser_getsecnkeys(dict, secName);
		for (k = 0; k < keysCount; ++k) {
//			printf("Checking %s...\n", keys[k]);
			found = 0;
			for (m = 0; m < expectedCount; ++m) {
				if (expected[m] != NULL && strcasecmp(keys[k], expected[m]) == 0) {
					found = 1;
					break;
				}
			}

			if (found == 0) {
				printf("WARNING: key %s from mcluster.ini is not recognized! Check the file mocca.ini for errors.\n", keys[k]);
			}
		}
	}

	iniparser_freedict(dict);
}


