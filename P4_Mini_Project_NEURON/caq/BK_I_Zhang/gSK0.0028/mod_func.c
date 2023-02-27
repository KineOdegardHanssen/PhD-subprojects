#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _CaDynamics_reg();
extern void _caq_reg();
extern void _kaf_reg();
extern void _naf_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," CaDynamics.mod");
fprintf(stderr," caq.mod");
fprintf(stderr," kaf.mod");
fprintf(stderr," naf.mod");
fprintf(stderr, "\n");
    }
_CaDynamics_reg();
_caq_reg();
_kaf_reg();
_naf_reg();
}
