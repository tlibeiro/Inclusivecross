#!/bin/csh
set tag="winter14"
cmsenv
rm AnalyzerIncHistosData276TeV_${tag}.root
rm AnalyzerIncHistosMC_${tag}.root
rm AnalyzerIncHistos_MC_${tag}_stabplots.root
./runAnalyzer
hadd AnalyzerIncHistosMC_${tag}.root AnalyzerIncHistos_allpt[0-9]_${tag}.root AnalyzerIncHistos_allpt10_${tag}.root
hadd AnalyzerIncHistos_MC_${tag}_stabplots.root AnalyzerIncHistos_allpt10_stabplots.root AnalyzerIncHistos_allpt[0-9]_stabplots.root
