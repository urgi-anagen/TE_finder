#/bin/bash

../NWalign seq1.fa seq2.fa > NWalign.outres
diff -s NWalign.outres NWalign.outexp
