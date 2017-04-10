all:

OPTFLAGS = -O0 -g -rdynamic -std=c++0x
OUTPUTDIR = ./

CFLAGS += $(OPTFLAGS) -fopenmp

#CXXFLAGS= -std=c++11 -D USEBYTE -D NEWHASH -D PREFETCH -g 

GITVERSION= "no-version"

CXXFLAGS+= -std=gnu++11 -DGITVERSION=\"$(GITVERSION)\"
CXXFLAGS+= -Wno-format -Wno-pointer-arith
CXXFLAGS+= -DKMERLENGTH=20   -march=native -DHASHSEED1=0x74b0dc51 -DHASHSEED2=0x19495cff -O3

CXXFLAGS+= -fopenmp

SRC= src/Ye_implementation/testOthelloKeyValue.cc

HEADERS= src/Ye_implementation/*.h

#all_test: test2 test4 test8 test16

assignMetagenomicsRead_allTaxoRank_12_w2: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D W2	

clean_assignMetagenomicsRead_allTaxoRank_12_w2:
	rm assignMetagenomicsRead_allTaxoRank_12_w2

assignMetagenomicsRead_allTaxoRank_12: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12	

clean_assignMetagenomicsRead_allTaxoRank_12:
	rm assignMetagenomicsRead_allTaxoRank_12


examineKaijuNodeDmpWithSimData:
	g++ $(CFLAGS) -o $(OUTPUTDIR)examineKaijuNodeDmpWithSimData src/metagenomics/utils/examineKaijuNodeDmpWithSimData.cpp

clean_examineKaijuNodeDmpWithSimData:
	rm examineKaijuNodeDmpWithSimData

reformKaijuResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformKaijuResults src/metagenomics/utils/reformKaijuResults.cpp

clean_reformKaijuResults:
	rm reformKaijuResults

reformKaijuResults_nodeDmp:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformKaijuResults_nodeDmp src/metagenomics/utils/reformKaijuResults_nodeDmp.cpp

clean_reformKaijuResults_nodeDmp:
	rm reformKaijuResults_nodeDmp	

parseTaxoNodeDmp2speciesId2taxoInfo:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseTaxoNodeDmp2speciesId2taxoInfo src/metagenomics/utils/parseTaxoNodeDmp2speciesId2taxoInfo.cpp

clean_parseTaxoNodeDmp2speciesId2taxoInfo:
	rm parseTaxoNodeDmp2speciesId2taxoInfo

assignRead2gene_16: src/geneCount/assign/assignRead2gene.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/geneCount/assign/assignRead2gene.cpp -DVALUELENGTH=16

assignRead2fusion_16: src/geneCount/assign/assignRead2gene.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/geneCount/assign/assignRead2gene.cpp -DVALUELENGTH=16 -D FUSION_GENE_READ_DETECTION

assignRead2gene_16_debug: src/geneCount/assign/assignRead2gene.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/geneCount/assign/assignRead2gene.cpp -DVALUELENGTH=16 -D PRINT_DEBUG

assignRead2fusion_16_debug: src/geneCount/assign/assignRead2gene.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/geneCount/assign/assignRead2gene.cpp -DVALUELENGTH=16 -D FUSION_GENE_READ_DETECTION -D PRINT_DEBUG

simulator:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulator src/geneCount/utils/simulator.cpp

countSimulatedGeneCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSimulatedGeneCount src/geneCount/utils/countSimulatedGeneCount.cpp	

getTaxoInfo_from_speciesFa2taxoInfoFile:
	g++ $(CXXFLAGS) -o $(OUTPUTDIR)getTaxoInfo_from_speciesFa2taxoInfoFile src/metagenomics/utils/getTaxoInfo_from_speciesFa2taxoInfoFile.cpp

binFloatList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)binFloatList src/utils/binFloatList.cpp

clean_binFloatList:
	rm binFloatList

countReadFileBaseNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countReadFileBaseNum src/utils/countReadFileBaseNum.cpp

clean_countReadFileBaseNum:
	rm countReadFileBaseNum

countReadFileBaseNum_withReadLengthMin:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countReadFileBaseNum_withReadLengthMin src/utils/countReadFileBaseNum_withReadLengthMin.cpp

clean_countReadFileBaseNum_withReadLengthMin:
	rm countReadFileBaseNum_withReadLengthMin

extractDiffAssignment:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractDiffAssignment src/metagenomics/utils/extractDiffAssignment.cpp

clean_extractDiffAssignment:
	rm extractDiffAssignment

summarizeTaxoAbundanceFromAssignmentResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarizeTaxoAbundanceFromAssignmentResults src/metagenomics/simulator/summarizeTaxoAbundanceFromAssignmentResults.cpp

clean_summarizeTaxoAbundanceFromAssignmentResults:
	rm summarizeTaxoAbundanceFromAssignmentResults

getSensiPreci_allTaxoRank:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSensiPreci_allTaxoRank src/metagenomics/simulator/getSensiPreci_allTaxoRank.cpp

clean_getSensiPreci_allTaxoRank:
	rm getSensiPreci_allTaxoRank

getSensiPreci_allTaxoRank_withReadLengthMin:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSensiPreci_allTaxoRank_withReadLengthMin src/metagenomics/simulator/getSensiPreci_allTaxoRank_withReadLengthMin.cpp

clean_getSensiPreci_allTaxoRank_withReadLengthMin:
	rm getSensiPreci_allTaxoRank_withReadLengthMin

getSensiPreci_allTaxoRank_withReadLengthMin_givenReadNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSensiPreci_allTaxoRank_withReadLengthMin_givenReadNum src/metagenomics/simulator/getSensiPreci_allTaxoRank_withReadLengthMin_givenReadNum.cpp

clean_getSensiPreci_allTaxoRank_withReadLengthMin_givenReadNum:
	rm getSensiPreci_allTaxoRank_withReadLengthMin_givenReadNum

combineAllTaxoRankAssignmentResults_Clark:
	g++ $(CFLAGS) -o $(OUTPUTDIR)combineAllTaxoRankAssignmentResults_Clark src/metagenomics/simulator/combineAllTaxoRankAssignmentResults_Clark.cpp

clean_combineAllTaxoRankAssignmentResults_Clark:
	rm combineAllTaxoRankAssignmentResults_Clark

getMappingRate_allTaxoRank:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMappingRate_allTaxoRank src/metagenomics/utils/getMappingRate_allTaxoRank.cpp	

clean_getMappingRate_allTaxoRank:
	rm getMappingRate_allTaxoRank

getMappingRate_specificTaxoRank_Clark:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMappingRate_specificTaxoRank_Clark src/metagenomics/utils/getMappingRate_specificTaxoRank_Clark.cpp

clean_getMappingRate_specificTaxoRank_Clark:
	rm getMappingRate_specificTaxoRank_Clark

merge_sortedMetagenomicsKmerSetIdFile_sortedHumanKmerCountFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)merge_sortedMetagenomicsKmerSetIdFile_sortedHumanKmerCountFile src/metagenomics/build/merge_sortedMetagenomicsKmerSetIdFile_sortedHumanKmerCountFile.cpp

clean_merge_sortedMetagenomicsKmerSetIdFile_sortedHumanKmerCountFile:
	rm merge_sortedMetagenomicsKmerSetIdFile_sortedHumanKmerCountFile

combine2AssignmentResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)combine2AssignmentResults src/metagenomics/assignRead/combine2AssignmentResults.cpp

clean_combine2AssignmentResults:
	rm combine2AssignmentResults

extractIncorrectRead_butOtherToolCorrect:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractIncorrectRead_butOtherToolCorrect src/metagenomics/utils/extractIncorrectRead_butOtherToolCorrect.cpp

clean_extractIncorrectRead_butOtherToolCorrect:
	rm extractIncorrectRead_butOtherToolCorrect

extractIncorrectRead:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractIncorrectRead src/metagenomics/utils/extractIncorrectRead.cpp

clean_extractIncorrectRead:
	rm extractIncorrectRead	

reformatKrakenSimBA5data:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatKrakenSimBA5data src/metagenomics/simulator/reformatKrakenSimBA5data.cpp

clean_reformatKrakenSimBA5data:
	rm reformatKrakenSimBA5data

reformatKrakenHiSeqMiSeqData:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatKrakenHiSeqMiSeqData src/metagenomics/simulator/reformatKrakenHiSeqMiSeqData.cpp

clean_reformatKrakenHiSeqMiSeqData:
	rm reformatKrakenHiSeqMiSeqData

evaluateAssignmentResults_Clark:
	g++ $(CFLAGS) -o $(OUTPUTDIR)evaluateAssignmentResults_Clark src/metagenomics/simulator/evaluateAssignmentResults_Clark.cpp

clean_evaluateAssignmentResults_Clark:
	rm evaluateAssignmentResults_Clark

evaluateAssignmentResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)evaluateAssignmentResults src/metagenomics/simulator/evaluateAssignmentResults.cpp

clean_evaluateAssignmentResults:
	rm evaluateAssignmentResults

evaluateAssignmentResults_allTaxoRankFormat:
	g++ $(CFLAGS) -o $(OUTPUTDIR)evaluateAssignmentResults_allTaxoRankFormat src/metagenomics/simulator/evaluateAssignmentResults_allTaxoRankFormat.cpp

clean_evaluateAssignmentResults_allTaxoRankFormat:
	rm evaluateAssignmentResults_allTaxoRankFormat

assignMetagenomicsRead_allTaxoRank_12_w2_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D W2	-D ASSIGN_INFO

clean_assignMetagenomicsRead_allTaxoRank_12_w2_assignInfo:
	rm assignMetagenomicsRead_allTaxoRank_12_w2_assignInfo


assignMetagenomicsRead_allTaxoRank_13_w2: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=13 -D W2	

clean_assignMetagenomicsRead_allTaxoRank_13_w2:
	rm assignMetagenomicsRead_allTaxoRank_13_w2

assignMetagenomicsRead_allTaxoRank_13_w2_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=13 -D W2	-D ASSIGN_INFO

clean_assignMetagenomicsRead_allTaxoRank_13_w2_assignInfo:
	rm assignMetagenomicsRead_allTaxoRank_13_w2_assignInfo


assignMetagenomicsRead_allTaxoRank_12_noFiltering: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D NO_FILTERING

clean_assignMetagenomicsRead_allTaxoRank_12_noFiltering:
	rm assignMetagenomicsRead_allTaxoRank_12_noFiltering

assignMetagenomicsRead_allTaxoRank_12_noFiltering_w2: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D NO_FILTERING -D NO_FILTERING_AT_BOTH_ENDS -D W2

clean_assignMetagenomicsRead_allTaxoRank_12_noFiltering_w2:
	rm assignMetagenomicsRead_allTaxoRank_12_noFiltering_w2



assignMetagenomicsRead_allTaxoRank_12_noFilteringAtBothEnds: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D NO_FILTERING_AT_BOTH_ENDS

clean_assignMetagenomicsRead_allTaxoRank_12_noFilteringAtBothEnds:
	rm assignMetagenomicsRead_allTaxoRank_12_noFilteringAtBothEnds	

assignMetagenomicsRead_allTaxoRank_12_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12	-D ASSIGN_INFO

clean_assignMetagenomicsRead_allTaxoRank_12_assignInfo:
	rm assignMetagenomicsRead_allTaxoRank_12_assignInfo

assignMetagenomicsRead_multiTaxoRank_12: src/metagenomics/assignRead/assignMetagenomicsRead_multiTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_multiTaxoRank.cpp -DVALUELENGTH=12	

clean_assignMetagenomicsRead_multiTaxoRank_12:
	rm assignMetagenomicsRead_multiTaxoRank_12

assignMetagenomicsRead_multiTaxoRank_12_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead_multiTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_multiTaxoRank.cpp -DVALUELENGTH=12	-D ASSIGN_INFO

clean_assignMetagenomicsRead_multiTaxoRank_12_assignInfo:
	rm assignMetagenomicsRead_multiTaxoRank_12_assignInfo

reformFqReadId:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformFqReadId src/metagenomics/simulator/reformFqReadId.cpp

clean_reformFqReadId:
	rm reformFqReadId

simulateReadForFaVec:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateReadForFaVec src/metagenomics/simulator/simulateReadForFaVec.cpp

clean_simulateReadForFaVec:
	rm simulateReadForFaVec

mergeSetSpecificKmerAtMultiTaxoLevel:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeSetSpecificKmerAtMultiTaxoLevel src/metagenomics/build/mergeSetSpecificKmerAtMultiTaxoLevel.cpp

clean_mergeSetSpecificKmerAtMultiTaxoLevel:
	rm mergeSetSpecificKmerAtMultiTaxoLevel

faVec2JfVec2KmerVec:
	g++ $(CFLAGS) -o $(OUTPUTDIR)faVec2JfVec2KmerVec src/utils/faVec2JfVec2KmerVec.cpp

clean_faVec2JfVec2KmerVec:
	rm faVec2JfVec2KmerVec

getSetSpecificKmerProfile_specificRank: src/metagenomics/getSetSpecificKmerProfile_specificRank.cpp
	g++ $(CXXFLAGS) -o $@ src/metagenomics/getSetSpecificKmerProfile_specificRank.cpp

clean_getSetSpecificKmerProfile_specificRank:
	rm getSetSpecificKmerProfile_specificRank

getTaxoKmerFreqProfile: src/metagenomics/getTaxoKmerFreqProfile.cpp
	g++ $(CXXFLAGS) -o $@ src/metagenomics/getTaxoKmerFreqProfile.cpp

clean_getTaxoKmerFreqProfile:
	rm getTaxoKmerFreqProfile

getTaxoSpecificKmerProfile: src/metagenomics/getTaxoSpecificKmerProfile.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ $< -DVALUELENGTH=16

rm_getTaxoSpecificKmerProfile:
	rm getTaxoSpecificKmerProfile

getSetSpecificKmerProfile: src/utils/getSetSpecificKmerProfile.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ $< -DVALUELENGTH=16

clean_getSetSpecificKmerProfile:
	rm getSetSpecificKmerProfile;

assignMetagenomicsRead_12: src/metagenomics/assignRead/assignMetagenomicsRead.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead.cpp -DVALUELENGTH=12	

clean_assignMetagenomicsRead_12:
	rm assignMetagenomicsRead_12;

assignMetagenomicsRead_12_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead.cpp -DVALUELENGTH=12 -D ASSIGN_INFO

clean_assignMetagenomicsRead_12_assignInfo:
	rm assignMetagenomicsRead_12_assignInfo

assignMetagenomicsRead_singleTaxoRank_12: src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp -DVALUELENGTH=12	

clean_assignMetagenomicsRead_singleTaxoRank_12:
	rm assignMetagenomicsRead_singleTaxoRank_12;

assignMetagenomicsRead_singleTaxoRank_12_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp -DVALUELENGTH=12 -D ASSIGN_INFO

clean_assignMetagenomicsRead_singleTaxoRank_12_assignInfo:
	rm assignMetagenomicsRead_singleTaxoRank_12_assignInfo

assignMetagenomicsRead_singleTaxoRank_6: src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp -DVALUELENGTH=6	

clean_assignMetagenomicsRead_singleTaxoRank_6:
	rm assignMetagenomicsRead_singleTaxoRank_6;

assignMetagenomicsRead_singleTaxoRank_6_assignInfo: src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignRead/assignMetagenomicsRead_singleTaxoRank.cpp -DVALUELENGTH=6 -D ASSIGN_INFO

clean_assignMetagenomicsRead_singleTaxoRank_6_assignInfo:
	rm assignMetagenomicsRead_singleTaxoRank_6_assignInfo		

test2: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=2

test4: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=4

test8: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=8

test16: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=16

clean_all_test:
	rm test2
	rm test4
	rm test8
	rm test16



buildKmerMPHindex: src/indexing/build/buildKmerMPHindex.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ $< -DVALUELENGTH=16

clean_buildKmerMPHindex:
	rm buildKmerMPHindex

KmerFileVec2KmerFrequencyProfile: src/utils/KmerFileVec2KmerFrequencyProfile.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ $< -DVALUELENGTH=16

clean_KmerFileVec2KmerFrequencyProfile:
	rm KmerFileVec2KmerFrequencyProfile

KmerFileVec2setSpecificAndRepetitiveKmerProfile: src/utils/KmerFileVec2setSpecificAndRepetitiveKmerProfile.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ $< -DVALUELENGTH=16

clean_KmerFileVec2setSpecificAndRepetitiveKmerProfile:
	rm KmerFileVec2setSpecificAndRepetitiveKmerProfile

KmerFileVec2setSpecificKmerProfile: src/indexing/building/KmerFileVec2setSpecificKmerProfile.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/indexing/building/KmerFileVec2setSpecificKmerProfile.cpp

#KmerFileVec2KmerFrequencyProfile: src/indexing/building/KmerFileVec2KmerFrequencyProfile.cpp src/Ye_implementation/*.h
#	g++ $(CXXFLAGS) -o $@ src/indexing/building/KmerFileVec2KmerFrequencyProfile.cpp	

#clean_KmerFileVec2setSpecificKmerProfile:
#	rm KmerFileVec2setSpecificKmerProfile

queryOthelloNode: queryOthelloNode_2 queryOthelloNode_4 queryOthelloNode_8 queryOthelloNode_16 

queryOthelloNode_2: src/queryOthelloNode.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/queryOthelloNode.cpp -DVALUELENGTH=2

queryOthelloNode_4: src/queryOthelloNode.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/queryOthelloNode.cpp -DVALUELENGTH=4

queryOthelloNode_8: src/queryOthelloNode.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/queryOthelloNode.cpp -DVALUELENGTH=8

queryOthelloNode_16: src/queryOthelloNode.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/queryOthelloNode.cpp -DVALUELENGTH=16

clean_all_queryOthelloNode:
	rm queryOthelloNode_2
	rm queryOthelloNode_4
	rm queryOthelloNode_8
	rm queryOthelloNode_16

buildLothelloTreeKmerProfileFromScratch:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildLothelloTreeKmerProfileFromScratch src/generateKmerProfile/buildLothelloTreeFromScratch.cpp

buildLothelloNodeKmerProfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildLothelloNodeKmerProfile src/generateKmerProfile/buildLothelloNodeKmerProfile.cpp

buildLothelloNodeBitArray_24: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=24

buildLothelloNodeBitArray_20: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=20

buildLothelloNodeBitArray_18: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=18

buildLothelloNodeBitArray_16: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=16

buildLothelloNodeBitArray_12: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=12	

buildLothelloNodeBitArray_8: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=8

buildLothelloNodeBitArray_6: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=6

buildLothelloNodeBitArray_4: src/generateBitArray/buildLothelloNodeBitArray.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/generateBitArray/buildLothelloNodeBitArray.cpp -DVALUELENGTH=4		

queryLothelloNode_Kmer_16: src/query/queryLothelloNode_Kmer.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_Kmer.cpp -DVALUELENGTH=16

queryLothelloNode_Kmer_12: src/query/queryLothelloNode_Kmer.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_Kmer.cpp -DVALUELENGTH=12

queryLothelloNode_Kmer_8: src/query/queryLothelloNode_Kmer.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_Kmer.cpp -DVALUELENGTH=8

queryLothelloNode_Kmer_6: src/query/queryLothelloNode_Kmer.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_Kmer.cpp -DVALUELENGTH=6

queryLothelloNode_Kmer_4: src/query/queryLothelloNode_Kmer.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_Kmer.cpp -DVALUELENGTH=4

queryLothelloNode_seq_8: src/query/queryLothelloNode_seq.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_seq.cpp -DVALUELENGTH=8

queryLothelloNode_seq_4: src/query/queryLothelloNode_seq.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_seq.cpp -DVALUELENGTH=4	

queryLothelloNode_file_8: src/query/queryLothelloNode_file.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_file.cpp -DVALUELENGTH=8

queryLothelloNode_file_4: src/query/queryLothelloNode_file.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/query/queryLothelloNode_file.cpp -DVALUELENGTH=4	

buildGenomeSubRegionSpecificKmerProfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildGenomeSubRegionSpecificKmerProfile src/generateKmerProfile/buildGenomeSubRegionSpecificKmerProfile.cpp

separateKmerCanonicalOrNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateKmerCanonicalOrNot src/generateKmerProfile/separateKmerCanonicalOrNot.cpp

getMergedGenome_fa:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMergedGenome_fa src/metagenomics/getMergedGenome_fa.cpp

getMergedSpecies_fa:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMergedSpecies_fa src/metagenomics/getMergedSpecies_fa.cpp	

getMergedGenome_jf:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMergedGenome_jf src/metagenomics/getMergedGenome_jf.cpp

getMergedGenome_jf_Kmer_sortedKmer:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMergedGenome_jf_Kmer_sortedKmer src/metagenomics/getMergedGenome_jf_Kmer_sortedKmer.cpp	

generateTaxoInfo_species_genus_phylum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateTaxoInfo_species_genus_phylum src/metagenomics/generateTaxoInfo_species_genus_phylum.cpp

readAssignment_4: src/genomeRegionAssignment/readAssignment.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/genomeRegionAssignment/readAssignment.cpp -DVALUELENGTH=4

readAssignment_6: src/genomeRegionAssignment/readAssignment.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/genomeRegionAssignment/readAssignment.cpp -DVALUELENGTH=6	

readAssignment_4_assignInfo: src/genomeRegionAssignment/readAssignment.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/genomeRegionAssignment/readAssignment.cpp -DVALUELENGTH=4 -D ASSIGN_INFO

readAssignment_6_assignInfo: src/genomeRegionAssignment/readAssignment.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/genomeRegionAssignment/readAssignment.cpp -DVALUELENGTH=6 -D ASSIGN_INFO	

extractGenomeRepeatRegion_4: src/genomeRegionAssignment/extractGenomeRepeatRegion.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/genomeRegionAssignment/extractGenomeRepeatRegion.cpp -DVALUELENGTH=4

mergeRepeatKmerWindow: src/genomeRegionAssignment/mergeRepeatKmerWindow.cpp src/Ye_implementation/*.h
	g++ $(CXXFLAGS) -o $@ src/genomeRegionAssignment/mergeRepeatKmerWindow.cpp -DVALUELENGTH=4	



KmerFile2sortedKmerFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)KmerFile2sortedKmerFile src/utils/KmerFile2sortedKmerFile.cpp

clean_KmerFile2sortedKmerFile:
	rm KmerFile2sortedKmerFile

KmerFile2sortedKmerFile_prefix3charOnly:
	g++ $(CFLAGS) -o $(OUTPUTDIR)KmerFile2sortedKmerFile_prefix3charOnly src/utils/KmerFile2sortedKmerFile_prefix3charOnly.cpp

clean_KmerFile2sortedKmerFile_prefix3charOnly:
	rm KmerFile2sortedKmerFile_prefix3charOnly

clean_buildLothelloTreeKmerProfileFromScratch:
	rm buildLothelloTreeKmerProfileFromScratch

clean_buildLothelloNodeKmerProfile:
	rm buildLothelloNodeKmerProfile

clean_buildLothelloNodeBitArray_8:
	rm LothelloNodeBitArray_8

clean_queryLothelloNode_Kmer_8:
	rm queryLothelloNode_singleKmer_8

clean_queryLothelloNode_seq_8:
	rm queryLothelloNode_seq_8

clean_buildGenomeSubRegionSpecificKmerProfile:
	rm buildGenomeSubRegionSpecificKmerProfile

clean_separateKmerCanonicalOrNot:
	rm separateKmerCanonicalOrNot

clean_getMergedGenome_fa:
	rm getMergedGenome_fa

clean_getMergedGenome_jf:
	rm getMergedGenome_jf	

clean_getMergedGenome_jf_Kmer_sortedKmer:
	rm getMergedGenome_jf_Kmer_sortedKmer

clean_generateTaxoInfo_species_genus_phylum:
	rm generateTaxoInfo_species_genus_phylum

clean_readAssignment_4:
	rm readAssignment_4

clean_readAssignment_4_assignInfo:
	rm readAssignment_4_assignInfo

clean_extractGenomeRepeatRegion_4:
	rm extractGenomeRepeatRegion_4

clean_simulator:
	rm simulator

clean_countSimulatedGeneCount:
	rm countSimulatedGeneCount

clean_assignRead2gene_16:
	rm assignRead2gene_16

clean_assignRead2fusion_16:
	rm assignRead2fusion_16	