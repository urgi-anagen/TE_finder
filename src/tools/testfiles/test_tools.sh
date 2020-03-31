#/bin/bash

function TestToolStdout
{
    TOOL=$1
    CMD=$2
    EXEC_DIR="../../../cmake-build-debug/src/tools"
    mkdir $TOOL
    cd $TOOL
    ../../$TOOL $EXEC_DIR/$CMD > $TOOL.outres
    if [ $? -ne 0 ]
    then
        echo "***PRG ERROR***"
    fi
    if [ ! -f ../$TOOL.outexp ]
    then
      echo $TOOL ": generate new expected output"
      cp $TOOL.outres ../$TOOL.outexp
    fi
    diff -q -s $TOOL.outres ../$TOOL.outexp
    if [ $? -ne 0 ]
    then
        echo "***OUTPUT ERROR***"
    fi
    cd ..
    rm -rf $TOOL
}

function TestTool
{
    TOOL=$1
    CMD=$2
    OUT=$3
    mkdir $TOOL
    cd $TOOL
    ../../$TOOL $CMD > $TOOL.outres
    if [ $? -ne 0 ]
    then
        echo "***PRG ERROR***"
    fi
    if [ ! -f ../$TOOL.outexp ]
    then
      echo $TOOL ": generate new expected output"
      cp $OUT ../$TOOL.outexp
    fi
    diff -q -s $OUT ../$TOOL.outexp
    if [ $? -ne 0 ]
    then
        echo "***OUTPUT ERROR***"
    fi
    cd ..
    rm -rf $TOOL
}
#--------------------

TestToolStdout NWalign "../seq1.fa ../seq2.fa"
TestToolStdout SWalign "../seq1.fa ../seq2.fa"
TestToolStdout TRsearch "../DmelChr4_refTEs.fa"
TestToolStdout polyAtail "../DmelChr4_refTEs.fa"
TestTool cutterDB "../DmelChr4.fa" ../DmelChr4.fa_cut
TestTool hsearch "-o out ../DmelChr4.fa ../DmelChr4_refTEs.fa" out
TestTool hrepeat "../DmelChr4_refTEs.fa" ../DmelChr4_refTEs.fa.hrepeat.set
TestToolStdout ltrsearch "../DmelChr4_refTEs.fa"
TestToolStdout fastlalign "../seq1.fa ../seq2.fa"
TestToolStdout rpt_map "../seq12.fa 10 -8 4 2"

#TestTool halign "../DmelChr4.fa ../DmelChr4_refTEs.fa"

#map2db      
#orienter    
#filterAlign 
#map2flank3  
#pathnum2id  
#halign      
#map2flank5  
#TSDsearch        
#map2flank53 
#refalign
#align2piler  
#mapOp       
#setnum2id
#align2recon 
#itrsearch
#mapdist     
#statalign   
#mapview.cp
#tabnum2id