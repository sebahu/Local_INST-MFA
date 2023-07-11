# Based upon Isosim, (C) 2019, INRA, France, GNU General Public License v3 (see license.txt for details)
# https://github.com/MetaSys-LISBP/IsoSim
#
# https://doi.org/10.1101/735308
# Path of the model

#from https://gist.github.com/jasonsychau/ff6bc78a33bf3fd1c6bd4fa78bbf42e7
stub <- function() {}
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
    }
}

setwd(thisPath())

setup_env <- function() {

    # Clean workspace
    rm(list= ls())
    
    library("vioplot")

    wd = getwd()
    # Load IsoSim
    setwd(file.path(dirname(dirname(wd)), "isosim"))
    source("isosim.R")

    # Go back to working directory and create "res" folder to store the results
    setwd(wd)
    if (!file.exists("res")){
        dir.create(file.path(wd, "res"))
    }
    setwd(file.path(wd, "res"))

    #########################
    ### GLOBAL PARAMETERS ###
    #########################
    
    # number of cores to use in parallel, i.e. at most how many child processes will be run simultaneously
    # single-core version can be used by setting 'numCores' to NULL
    .GlobalEnv$numCores <- detectCores()
    
    # number of Monte Carlo iterations for flux calculation
    .GlobalEnv$mc_iter <- 100
    
}

init_aracore_model <- function() {
    
    ####################################
    cat("\n   ... construct isotopic model of AraCore ...\n\n")
    ####################################
    
    # network definition
    # reactions: Reactants, Coefficients, Rate law, Atom transitions
    .GlobalEnv$rxn <- list(r5     = list("R"=c("ADPzzh" ,"ATPzzh"), "C"=c(-3 ,3),     "E"="v5",   "T"=c("ABCDE" ,"ABCDE")),
                r7     = list("R"=c("ADPzzh" ,"ATPzzh"), "C"=c(1 ,-1),     "E"="v7",   "T"=c("ABCDE" ,"ABCDE")),
                r18     = list("R"=c("ADPzzh" ,"ATPzzh"), "C"=c(1 ,-1),     "E"="v18",   "T"=c("ABCDE" ,"ABCDE")),
                r21     = list("R"=c("ATPzzh" ,"ADPGzzh"), "C"=c(-1 ,1),     "E"="v21",   "T"=c("ABCDE" ,"ABCDE")),
                r22     = list("R"=c("ADPzzh" ,"ADPGzzh"), "C"=c(1 ,-1),     "E"="v22",   "T"=c("ABCDE" ,"ABCDE")),
                r23     = list("R"=c("ADPzzh" ,"ADPGzzh"), "C"=c(1 ,-1),     "E"="v23",   "T"=c("ABCDE" ,"ABCDE")),
                r24     = list("R"=c("ADPzzh" ,"ADPGzzh"), "C"=c(1 ,-1),     "E"="v24",   "T"=c("ABCDE" ,"ABCDE")),
                r41     = list("R"=c("UTPzzc" ,"UDPGzzc"), "C"=c(-1 ,1),     "E"="v41",   "T"=c("AB" ,"AB")),
                r42     = list("R"=c("UDPGzzc" ,"UDPzzc"), "C"=c(-1 ,1),     "E"="v42",   "T"=c("AB" ,"AB")),
                r44     = list("R"=c("UDPGzzc" ,"UDPzzc"), "C"=c(-1 ,1),     "E"="v44",   "T"=c("AB" ,"AB")),
                r47     = list("R"=c("UDPGzzc" ,"UDPzzc"), "C"=c(-1 ,1),     "E"="v47",   "T"=c("AB" ,"AB")),
                r49     = list("R"=c("UDPGzzc" ,"UDPzzc"), "C"=c(-1 ,1),     "E"="v49",   "T"=c("AB" ,"AB")),
                r50     = list("R"=c("UDPGzzc" ,"UDPzzc"), "C"=c(-1 ,1),     "E"="v50",   "T"=c("AB" ,"AB")),
                r52     = list("R"=c("NADzzc" ,"NADHzzc"), "C"=c(1 ,-1),     "E"="v52",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r57     = list("R"=c("ATPzzc" ,"ADPzzc"), "C"=c(1 ,-1),     "E"="v57",   "T"=c("ABCDE" ,"ABCDE")),
                r58     = list("R"=c("NADzzc" ,"NADHzzc"), "C"=c(-1 ,1),     "E"="v58",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r59     = list("R"=c("NADPzzc" ,"NADPHzzc"), "C"=c(-1 ,1),     "E"="v59",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r62     = list("R"=c("NADzzh" ,"NADHzzh"), "C"=c(1 ,-1),     "E"="v62",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r67     = list("R"=c("ThPPzzm" ,"HxxEthxxThPPzzm"), "C"=c(-1 ,1),     "E"="v67",   "T"=c("ABCD" ,"ABCD")),
                r68     = list("R"=c("ThPPzzm" ,"HxxEthxxThPPzzm" ,"LPAzzm" ,"AxxDHLzzm"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v68",   "T"=c("ABCD" ,"ABCD" ,"E" ,"E")),
                r69     = list("R"=c("AxxDHLzzm" ,"CoAzzm" ,"AxxCoAzzm" ,"DHLzzm"), "C"=c(-1 ,-1 ,1 ,1),     "E"="v69",   "T"=c("A" ,"BCDEFGH" ,"BCDEFGH" ,"A")),
                r70     = list("R"=c("LPAzzm" ,"DHLzzm" ,"NADzzm" ,"NADHzzm"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v70",   "T"=c("A" ,"A" ,"BCDEFGH" ,"BCDEFGH")),
                r71     = list("R"=c("CoAzzm" ,"AxxCoAzzm"), "C"=c(1 ,-1),     "E"="v71",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r75     = list("R"=c("LPAzzm" ,"SxxDHLzzm"), "C"=c(-1 ,1),     "E"="v75",   "T"=c("A" ,"A")),
                r76     = list("R"=c("CoAzzm" ,"DHLzzm" ,"SxxDHLzzm" ,"SxxCoAzzm"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v76",   "T"=c("ABCDEFG" ,"H" ,"H" ,"ABCDEFG")),
                r77     = list("R"=c("CoAzzm" ,"SxxCoAzzm" ,"ADPzzm" ,"ATPzzm"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v77",   "T"=c("FGHIJKL" ,"FGHIJKL" ,"ABCDE" ,"ABCDE")),
                r80     = list("R"=c("NADzzm" ,"NADHzzm"), "C"=c(-1 ,1),     "E"="v80",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r81     = list("R"=c("NADzzm" ,"NADHzzm"), "C"=c(1 ,-1),     "E"="v81",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r84     = list("R"=c("ADPzzm" ,"ATPzzm"), "C"=c(-1 ,1),     "E"="v84",   "T"=c("ABCDE" ,"ABCDE")),
                r90     = list("R"=c("Glyzzm" ,"LPLzzm" ,"amDHPzzm"), "C"=c(-1 ,-1 ,1),     "E"="v90",   "T"=c("A" ,"B" ,"AB")),
                r91     = list("R"=c("amDHPzzm" ,"THFzzm" ,"DHPzzm" ,"MxxTHFzzm" ,"NH4zzm"), "C"=c(-1 ,-1 ,1 ,1 ,1),     "E"="v91",   "T"=c("AB" ,"CDEFGHI" ,"B" ,"CEGHIDF" ,"A")),
                r92     = list("R"=c("NADzzm" ,"NADHzzm" ,"LPLzzm" ,"DHPzzm"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v92",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r93     = list("R"=c("Glyzzm" ,"THFzzm" ,"MxxTHFzzm" ,"Serzzm"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v93",   "T"=c("A" ,"BGCHDEF" ,"BCDEFGH" ,"A")),
                r100     = list("R"=c("NADzzh" ,"NADHzzh"), "C"=c(-1 ,1),     "E"="v100",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r101     = list("R"=c("NADPzzh" ,"NADPHzzh"), "C"=c(1 ,-1),     "E"="v101",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r106     = list("R"=c("NADPzzh" ,"NADPHzzh"), "C"=c(1 ,-1),     "E"="v106",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r107     = list("R"=c("ADPzzh" ,"ATPzzh"), "C"=c(1 ,-1),     "E"="v107",   "T"=c("ABCDE" ,"ABCDE")),
                r110     = list("R"=c("ATPzzh" ,"AMPzzh"), "C"=c(-1 ,1),     "E"="v110",   "T"=c("ABCDE" ,"ABCDE")),
                r111     = list("R"=c("ATPzzh" ,"AMPzzh"), "C"=c(1 ,-1),     "E"="v111",   "T"=c("ABCDE" ,"ABCDE")),
                r113     = list("R"=c("NADPzzh" ,"NADPHzzh"), "C"=c(-1 ,1),     "E"="v113",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r114     = list("R"=c("NADzzh" ,"NADHzzh"), "C"=c(-1 ,1),     "E"="v114",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r124     = list("R"=c("NADzzc" ,"NADHzzc"), "C"=c(-1 ,1),     "E"="v124",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r129     = list("R"=c("AxxCoAzzc" ,"CoAzzc"), "C"=c(-1 ,1),     "E"="v129",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r130     = list("R"=c("ATPzzc" ,"ADPzzc" ,"AxxCoAzzc" ,"CoAzzc"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v130",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJKL" ,"FGHIJKL")),
                r131     = list("R"=c("ATPzzh" ,"AMPzzh" ,"AxxCoAzzh" ,"CoAzzh"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v131",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJKL" ,"FGHIJKL")),
                r132     = list("R"=c("ADPzzh" ,"ATPzzh" ,"AxxCoAzzh" ,"MxxCoAzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v132",   "T"=c("HIJKL" ,"HIJKL" ,"ABCDEFG" ,"ABCDEFG")),
                r133     = list("R"=c("CoAzzh" ,"MxxCoAzzh"), "C"=c(1 ,-1),     "E"="v133",   "T"=c("CDEFGHI" ,"CDEFGHI" )),
                r137     = list("R"=c("NADzzh" ,"NADHzzh"), "C"=c(1 ,-1),     "E"="v137",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r138     = list("R"=c("ADPzzh" ,"ATPzzh" ,"THFzzh" ,"FxxTHFzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v138",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJKL" ,"FHIJKLG")),
                r139     = list("R"=c("NADPzzc" ,"NADPHzzc" ,"DHFzzc" ,"THFzzc"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v139",   "T"=c("HIJKLMN" ,"HIJKLMN" ,"ABCDEFG" ,"ABCDEFG")),
                r140     = list("R"=c("NADPzzc" ,"NADPHzzc"), "C"=c(1 ,-1),     "E"="v140",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r148     = list("R"=c("ATPzzh" ,"APSzzh"), "C"=c(-1 ,1),     "E"="v148",   "T"=c("ABCDE" ,"ABCDE")),
                r149     = list("R"=c("AMPzzh" ,"APSzzh" ), "C"=c(1 ,-1 ),     "E"="v149",   "T"=c("ABCDE" ,"ABCDE" )),
                r150     = list("R"=c("NADPzzh" ,"NADPHzzh" ), "C"=c(1 ,-1 ),     "E"="v150",   "T"=c("ABCDEFG" ,"ABCDEFG" )),
                r155     = list("R"=c("Gluzzm" ,"Alazzm"), "C"=c(-1 ,1),     "E"="v155",   "T"=c("A" ,"A")),
                r156     = list("R"=c("Gluzzp" ,"Alazzp"), "C"=c(-1 ,1),     "E"="v156",   "T"=c("A" ,"A")),
                r157     = list("R"=c("ADPzzh" ,"ATPzzh" ,"AxxGluzzh" ,"AxxGluPzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v157",   "T"=c("BCDEF" ,"BCDEF" ,"A" ,"A")),
                r158     = list("R"=c("NADPzzh" ,"NADPHzzh" ,"AxxGluPzzh" ,"AxxGluxxSeAzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v158",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r159     = list("R"=c("Gluzzh" ,"AxxGluxxSeAzzh" ,"AxxOrnzzh"), "C"=c(-1 ,-1 ,1),     "E"="v159",   "T"=c("B" ,"A" ,"BA")),
                r160     = list("R"=c("Gluzzh" ,"AxxGluzzh" ,"AxxOrnzzh" ,"Ornzzh"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v160",   "T"=c("C" ,"C" ,"AB" ,"AB")),
                r161     = list("R"=c("ADPzzh" ,"ATPzzh" ,"Gluzzh" ,"Glnzzh" ,"CBPzzh"), "C"=c(2 ,-2 ,1 ,-1 ,1),     "E"="v161",   "T"=c("ABCDE" ,"ABCDE" ,"F" ,"FG" ,"G")),
                r162     = list("R"=c("Ornzzh" ,"CBPzzh" ,"CTLzzh"), "C"=c(-1 ,-1 ,1),     "E"="v162",   "T"=c("BC" ,"A" ,"CAB")),
                r163     = list("R"=c("ATPzzh" ,"AMPzzh" ,"CTLzzh" ,"Aspzzh" ,"ArgxxSCAzzh"), "C"=c(-1 ,1 ,-1 ,-1 ,1),     "E"="v163",   "T"=c("BCDEF" ,"BCDEF" ,"GHI" ,"A" ,"GHIA")),
                r164     = list("R"=c("ArgxxSCAzzh" ,"Argzzh"), "C"=c(-1 ,1),     "E"="v164",   "T"=c("ABCD" ,"ADBC")),
                r166     = list("R"=c("ATPzzc" ,"AMPzzc" ,"Glnzzc" ,"Aspzzc" ,"Gluzzc" ,"Asnzzc"), "C"=c(-1 ,1 ,-1 ,-1 ,1 ,1),     "E"="v166",   "T"=c("BCDEF" ,"BCDEF" ,"GH" ,"A" ,"G" ,"AH")),
                r167     = list("R"=c("Aspzzc" ,"Gluzzc"), "C"=c(1 ,-1),     "E"="v167",   "T"=c("A" ,"A")),
                r168     = list("R"=c("Gluzzh" ,"Aspzzh"), "C"=c(-1 ,1),     "E"="v168",   "T"=c("A" ,"A")),
                r172     = list("R"=c("AxxCoAzzc" ,"CoAzzc" ,"Serzzc" ,"AxxSerzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v172",   "T"=c("ABCDEFG" ,"ABCDEFG" ,"H" ,"H")),
                r173     = list("R"=c("AxxCoAzzh" ,"CoAzzh" ,"Serzzh" ,"AxxSerzzh"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v173",   "T"=c("ABCDEFG" ,"ABCDEFG" ,"H" ,"H")),
                r175     = list("R"=c("AxxSerzzc" ,"Cyszzc"), "C"=c(-1 ,1),     "E"="v175",   "T"=c("A" ,"A")),
                r176     = list("R"=c("AxxSerzzh" ,"Cyszzh"), "C"=c(-1 ,1),     "E"="v176",   "T"=c("A" ,"A")),
                r178     = list("R"=c("NADPzzh" ,"NADPHzzh" ,"Gluzzh" ,"GluBzzh" ,"Glnzzh"), "C"=c(-1 ,1 ,-1 ,-1 ,1),     "E"="v178",   "T"=c("CDEFGHI" ,"CDEFGHI" ,"A" , "B", "AB")),
                r178b   = list("R"=c("Gluzzh" ,"GluBzzh"), "C"=c(-1, 1),     "E"="v178b",   "T"=c("A" , "A")),
                r179     = list("R"=c("NADPzzc" ,"NADPHzzc" ,"Gluzzc" ,"NH4zzc"), "C"=c(1 ,-1 ,1 ,-1),     "E"="v179",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r182     = list("R"=c("Gluzzh" ,"Glnzzh"), "C"=c(2 ,-1),     "E"="v182",   "T"=c("A" ,"AB")),
                r184     = list("R"=c("NH4zzm" ,"NADPHzzm" ,"NADPzzm" ,"Gluzzm"), "C"=c(1 ,1 ,-1 ,-1),     "E"="v184",   "T"=c("H" ,"ABCDEFG" ,"ABCDEFG" ,"H")),
                r185     = list("R"=c("NADzzm" ,"NADHzzm" ,"NH4zzm" ,"Gluzzm"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v185",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r186     = list("R"=c("Gluzzc" ,"GABAzzc"), "C"=c(-1 ,1),     "E"="v186",   "T"=c("A" ,"A")),
                r187     = list("R"=c("Alazzm" ,"GABAzzm"), "C"=c(-1 ,1),     "E"="v187",   "T"=c("A" ,"A")),
                r188     = list("R"=c("Gluzzm" ,"GABAzzm"), "C"=c(1 ,-1),     "E"="v188",   "T"=c("A" ,"A")),
                r189     = list("R"=c("NADzzm" ,"NADHzzm"), "C"=c(1 ,-1),     "E"="v189",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r190     = list("R"=c("ATPzzc" ,"ADPzzc" ,"Glnzzc" ,"Gluzzc" ,"NH4zzc"), "C"=c(-1 ,1 ,1 ,-1 ,-1),     "E"="v190",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"F" ,"G")),
                r191     = list("R"=c("ADPzzh" ,"ATPzzh" ,"NH4zzh" ,"Gluzzh" ,"Glnzzh"), "C"=c(1 ,-1 ,-1 ,-1 ,1),     "E"="v191",   "T"=c("ABCDE" ,"ABCDE" ,"G" ,"F" ,"FG")),
                r193     = list("R"=c("Thrzzc" ,"Glyzzc"), "C"=c(-1 ,1),     "E"="v193",   "T"=c("A" ,"A")),
                r194     = list("R"=c("NADzzc" ,"NADHzzc"), "C"=c(-1 ,1),     "E"="v194",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r199     = list("R"=c("Gluzzh" ,"Glyzzh"), "C"=c(-1 ,1),     "E"="v199",   "T"=c("A" ,"A")),
                r200     = list("R"=c("ATPzzh" ,"PRxxATPzzh"), "C"=c(-1 ,1),     "E"="v200",   "T"=c("ABCDE" ,"ADCBE")),
                r201     = list("R"=c("PRxxATPzzh" ,"PRxxAMPzzh"), "C"=c(-1 ,1),     "E"="v201",   "T"=c("ABCDE" ,"ABCDE")),
                r202     = list("R"=c("PRxxAMPzzh" ,"PxxAICARxxPzzh"), "C"=c(-1 ,1),     "E"="v202",   "T"=c("ABCDE" ,"ACDBE")),
                r203     = list("R"=c("PxxAICARxxPzzh" ,"PuxxAICARxxPzzh"), "C"=c(-1 ,1),     "E"="v203",   "T"=c("ABCDE" ,"ACBDE")),
                r204     = list("R"=c("Gluzzh" ,"Glnzzh" ,"PuxxAICARxxPzzh" ,"EIGPzzh" ,"AICARzzh"), "C"=c(1 ,-1 ,-1 ,1 ,1),     "E"="v204",   "T"=c("A" ,"AB" ,"CDEFG" ,"DB" ,"ECFG")),
                r205     = list("R"=c("EIGPzzh" ,"IAxxPzzh"), "C"=c(-1 ,1),     "E"="v205",   "T"=c("AB" ,"AB")),
                r206     = list("R"=c("Gluzzh" ,"IAxxPzzh" ,"HisolxxPzzh"), "C"=c(-1 ,-1 ,1),     "E"="v206",   "T"=c("A" ,"BC" ,"ABC")),
                r207     = list("R"=c("HisolxxPzzh" ,"Hisolzzh"), "C"=c(-1 ,1),     "E"="v207",   "T"=c("ABC" ,"ABC")),
                r208     = list("R"=c("NADzzh" ,"NADHzzh" ,"Hisolzzh" ,"Hisalzzh"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v208",   "T"=c("DEFGHIJ" ,"DEFGHIJ" ,"ABC" ,"ABC")),
                r209     = list("R"=c("NADzzh" ,"NADHzzh" ,"Hisalzzh" ,"Hiszzh"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v209",   "T"=c("DEFGHIJ" ,"DEFGHIJ" ,"ABC" ,"ABC")),
                r210     = list("R"=c("NH4zzh" ,"Thrzzh"), "C"=c(1 ,-1),     "E"="v210",   "T"=c("A" ,"A")),
                r211     = list("R"=c("ThPPzzh" ,"HxxEthxxThPPzzh"), "C"=c(-1 ,1),     "E"="v211",   "T"=c("ABCD" ,"ABCD")),
                r212     = list("R"=c("ThPPzzh" ,"HxxEthxxThPPzzh"), "C"=c(1 ,-1),     "E"="v212",   "T"=c("ABCD" ,"ABCD")),
                r214     = list("R"=c("NADPzzh" ,"NADPHzzh"), "C"=c(1 ,-1),     "E"="v214",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r216     = list("R"=c("Gluzzh" ,"Ilezzh"), "C"=c(-1 ,1),     "E"="v216",   "T"=c("A" ,"A")),
                r217     = list("R"=c("ThPPzzh" ,"HxxEthxxThPPzzh"), "C"=c(1 ,-1),     "E"="v217",   "T"=c("ABCD" ,"ABCD")),
                r219     = list("R"=c("NADPzzh" ,"NADPHzzh"), "C"=c(1 ,-1),     "E"="v219",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r221     = list("R"=c("AxxCoAzzh" ,"CoAzzh"), "C"=c(-1 ,1),     "E"="v221",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r224     = list("R"=c("NADzzh" ,"NADHzzh"), "C"=c(-1 ,1),     "E"="v224",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r226     = list("R"=c("Gluzzh" ,"Leuzzh"), "C"=c(-1 ,1),     "E"="v226",   "T"=c("A" ,"A")),
                r227     = list("R"=c("ADPzzh" ,"ATPzzh" ,"Aspzzh" ,"AspPzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v227",   "T"=c("BCDEF" ,"BCDEF" ,"A" ,"A")),
                r228     = list("R"=c("NADPzzh" ,"NADPHzzh" ,"AspPzzh" ,"AspxxSeAzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v228",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r229     = list("R"=c("AspxxSeAzzh" ,"DHDzzh"), "C"=c(-1 ,1),     "E"="v229",   "T"=c("A" ,"A")),
                r231     = list("R"=c("NADPzzh" ,"NADPHzzh" ,"DHDzzh" ,"THDzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v231",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r232     = list("R"=c("Gluzzh" ,"THDzzh" ,"DAPzzh"), "C"=c(-1 ,-1 ,1),     "E"="v232",   "T"=c("A" ,"B" ,"AB")),
                r233     = list("R"=c("DAPzzh" ,"mDAPzzh"), "C"=c(-1 ,1),     "E"="v233",   "T"=c("AB" ,"AB")),
                r234     = list("R"=c("mDAPzzh" ,"Lyszzh"), "C"=c(-1 ,1),     "E"="v234",   "T"=c("AB" ,"BA")),
                r235     = list("R"=c("NADPzzh" ,"NADPHzzh" ,"AspxxSeAzzh" ,"HxxSerzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v235",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r237     = list("R"=c("ADPzzh" ,"ATPzzh" ,"HxxSerzzh" ,"PHxxSerzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v237",   "T"=c("ABCDE" ,"ABCDE" ,"F" ,"F")),
                r238     = list("R"=c("Cyszzh" ,"PHxxSerzzh" ,"CTHzzh"), "C"=c(-1 ,-1 ,1),     "E"="v238",   "T"=c("A" ,"B" ,"BA")),
                r239     = list("R"=c("NH4zzh" ,"CTHzzh" ,"HxxCyszzh"), "C"=c(1 ,-1 ,1),     "E"="v239",   "T"=c("B" ,"AB" ,"A")),
                r240     = list("R"=c("NADzzc" ,"NADHzzc" ,"MxxTHFzzc" ,"5MxxTHFzzc"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v240",   "T"=c("HIJKLMN" ,"HIJKLMN" ,"ABCDEFG" ,"AFBCDEG")),
                r241     = list("R"=c("THFzzc" ,"5MxxTHFzzc" ,"HxxCyszzc" ,"Metzzc"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v241",   "T"=c("ABCGDEF" ,"ABCDEFG" ,"H" ,"H")),
                r242     = list("R"=c("THFzzc" ,"Serzzc" ,"Glyzzc" ,"MxxTHFzzc"), "C"=c(-1 ,-1 ,1 ,1),     "E"="v242",   "T"=c("BGCHDEF" ,"A" ,"A" ,"BCDEFGH")),
                r247     = list("R"=c("Gluzzh" ,"AGNzzh"), "C"=c(-1 ,1),     "E"="v247",   "T"=c("A" ,"A")),
                r248     = list("R"=c("AGNzzh" ,"Phezzh"), "C"=c(-1 ,1),     "E"="v248",   "T"=c("A" ,"A")),
                r249     = list("R"=c("Argzzm" ,"Ornzzm" ,"ureazzm"), "C"=c(-1 ,1 ,1),     "E"="v249",   "T"=c("ABCD" ,"DA" ,"BC")),
                r251     = list("R"=c("Gluzzm" ,"Ornzzm" ,"GluxxSeAzzm"), "C"=c(1 ,-1 ,1),     "E"="v251",   "T"=c("A" ,"AB" ,"B")),
                r252     = list("R"=c("ATPzzc" ,"ADPzzc" ,"Gluzzc" ,"GluPzzc"), "C"=c(1 ,-1 ,1 ,-1),     "E"="v252",   "T"=c("ABCDE" ,"ABCDE" ,"F" ,"F")),
                r254     = list("R"=c("ADPzzm" ,"ATPzzm" ,"Gluzzm" ,"GluPzzm"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v254",   "T"=c("ABCDE" ,"ABCDE" ,"F" ,"F")),
                r255     = list("R"=c("NADPzzc" ,"NADPHzzc" ,"GluPzzc" ,"GluxxSeAzzc"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v255",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r257     = list("R"=c("NADPHzzm" ,"NADPzzm" ,"GluxxSeAzzm" ,"GluPzzm"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v257",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r258     = list("R"=c("GluxxSeAzzc" ,"P5Czzc"), "C"=c(1 ,-1),     "E"="v258",   "T"=c("A" ,"A")),
                r260     = list("R"=c("GluxxSeAzzm" ,"P5Czzm"), "C"=c(-1 ,1),     "E"="v260",   "T"=c("A" ,"A")),
                r261     = list("R"=c("NADzzc" ,"NADHzzc" ,"P5Czzc" ,"Prozzc"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v261",   "T"=c("ABCDEFG" ,"ABCDEFG" ,"H" ,"H")),
                r262     = list("R"=c("NADzzm" ,"NADHzzm" ,"P5Czzm" ,"Prozzm"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v262",   "T"=c("ABCDEFG" ,"ABCDEFG" ,"H" ,"H")),
                r267     = list("R"=c("NADzzh" ,"NADHzzh"), "C"=c(-1 ,1),     "E"="v267",   "T"=c("ABCDEFG" ,"ABCDEFG")),
                r268     = list("R"=c("Gluzzh" ,"PSerzzh"), "C"=c(-1 ,1),     "E"="v268",   "T"=c("A" ,"A")),
                r269     = list("R"=c("Serzzh" ,"PSerzzh"), "C"=c(1 ,-1),     "E"="v269",   "T"=c("A" ,"A")),
                r270     = list("R"=c("NH4zzc" ,"Serzzc"), "C"=c(1 ,-1),     "E"="v270",   "T"=c("A" ,"A")),
                r271     = list("R"=c("Thrzzh" ,"PHxxSerzzh"), "C"=c(1 ,-1),     "E"="v271",   "T"=c("A" ,"A")),
                r272     = list("R"=c("Gluzzh" ,"Glnzzh" ,"ANTzzh"), "C"=c(1 ,-1 ,1),     "E"="v272",   "T"=c("A" ,"AB" ,"B")),
                r273     = list("R"=c("ANTzzh" ,"PRxxANTzzh"), "C"=c(-1 ,1),     "E"="v273",   "T"=c("A" ,"A")),
                r274     = list("R"=c("PRxxANTzzh" ,"CPDxxRu5Pzzh"), "C"=c(-1 ,1),     "E"="v274",   "T"=c("A" ,"A")),
                r275     = list("R"=c("CPDxxRu5Pzzh" ,"IndxxGPzzh"), "C"=c(-1 ,1),     "E"="v275",   "T"=c("A" ,"A")),
                r276     = list("R"=c("IndxxGPzzh" ,"Indzzh"), "C"=c(-1 ,1),     "E"="v276",   "T"=c("A" ,"A")),
                r277     = list("R"=c("Serzzh" ,"Indzzh" ,"Trpzzh"), "C"=c(-1 ,-1 ,1),     "E"="v277",   "T"=c("B" ,"A" ,"BA")),
                r278     = list("R"=c("NADPzzh" ,"NADPHzzh" ,"AGNzzh" ,"Tyrzzh"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v278",   "T"=c("BCDEFGH" ,"BCDEFGH" ,"A" ,"A")),
                r279     = list("R"=c("Gluzzh" ,"Valzzh"), "C"=c(-1 ,1),     "E"="v279",   "T"=c("A" ,"A")),
                r280     = list("R"=c("Gluzzh" ,"Glnzzh" ,"PRAzzh"), "C"=c(1 ,-1 ,1),     "E"="v280",   "T"=c("A" ,"AB" ,"B")),
                r281     = list("R"=c("ADPzzh" ,"ATPzzh" ,"Glyzzh" ,"PRAzzh" ,"GARzzh"), "C"=c(1 ,-1 ,-1 ,-1 ,1),     "E"="v281",   "T"=c("ABCDE" ,"ABCDE" ,"F" ,"G" ,"FG")),
                r282     = list("R"=c("THFzzh" ,"FxxTHFzzh" ,"GARzzh" ,"FGARzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v282",   "T"=c("AGBCDEF" ,"ABCDEFG" ,"HI" ,"HI")),
                r283     = list("R"=c("ADPzzh" ,"ATPzzh" ,"Gluzzh" ,"Glnzzh" ,"FGARzzh" ,"FGAMzzh"), "C"=c(1 ,-1 ,1 ,-1 ,-1 ,1),     "E"="v283",   "T"=c("ABCDE" ,"ABCDE" ,"H" ,"HI" ,"FG" ,"IFG")),
                r284     = list("R"=c("ADPzzh" ,"ATPzzh" ,"FGAMzzh" ,"AIRzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v284",   "T"=c("ABCDE" ,"ABCDE" ,"FGH" ,"FGH")),
                r285     = list("R"=c("AIRzzh" ,"CAIRzzh"), "C"=c(-1 ,1),     "E"="v285",   "T"=c("ABC" ,"ABC")),
                r286     = list("R"=c("ADPzzh" ,"ATPzzh" ,"Aspzzh" ,"CAIRzzh" ,"SAICARzzh"), "C"=c(1 ,-1 ,-1 ,-1 ,1),     "E"="v286",   "T"=c("BCDEF" ,"BCDEF" ,"A" ,"GHI" ,"GHAI")),
                r287     = list("R"=c("AICARzzh" ,"SAICARzzh"), "C"=c(1 ,-1),     "E"="v287",   "T"=c("ACBD" ,"ABCD")),
                r288     = list("R"=c("THFzzh" ,"FxxTHFzzh" ,"AICARzzh" ,"FAICARzzh"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v288",   "T"=c("EKFGHIJ" ,"EFGHIJK" ,"ABCD" ,"BCAD")),
                r289     = list("R"=c("FAICARzzh" ,"IMPzzh"), "C"=c(-1 ,1),     "E"="v289",   "T"=c("ABCD" ,"CABD")),
                r290     = list("R"=c("AMPzzc" ,"NH4zzc" ,"IMPzzc"), "C"=c(-1 ,1 ,1),     "E"="v290",   "T"=c("ABCDE" ,"A" ,"CBDE")),
                r291     = list("R"=c("NADzzc" ,"NADHzzc" ,"IMPzzc" ,"XMPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v291",   "T"=c("EFGHIJK" ,"EFGHIJK" ,"ABCD" ,"CABD")),
                r292     = list("R"=c("ATPzzc" ,"AMPzzc" ,"Glnzzc" ,"Gluzzc" ,"XMPzzc" ,"GMPzzc"), "C"=c(-1 ,1 ,-1 ,1 ,-1 ,1),     "E"="v292",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"F" ,"HIJK" ,"GHIJK")),
                r293     = list("R"=c("Aspzzh" ,"IMPzzh" ,"GTPzzh" ,"DCxxAMPzzh" ,"GDPzzh"), "C"=c(-1 ,-1 ,-1 ,1 ,1),     "E"="v293",   "T"=c("A" ,"GHIJ" ,"BCDEF" ,"GHIAJ" ,"BCDEF")),
                r294     = list("R"=c("AMPzzh" ,"DCxxAMPzzh"), "C"=c(1 ,-1),     "E"="v294",   "T"=c("DABCE" ,"ABCDE")),
                r295     = list("R"=c("CBPzzh" ,"Aspzzh" ,"CBxxAspzzh"), "C"=c(-1 ,-1 ,1),     "E"="v295",   "T"=c("B" ,"A" ,"BA")),
                r296     = list("R"=c("CBxxAspzzh" ,"DHOzzh"), "C"=c(-1 ,1),     "E"="v296",   "T"=c("AB" ,"BA")),
                r297     = list("R"=c("DHOzzm" ,"OROzzm"), "C"=c(-1 ,1),     "E"="v297",   "T"=c("AB" ,"AB")),
                r298     = list("R"=c("OROzzc" ,"OMPzzc"), "C"=c(-1 ,1),     "E"="v298",   "T"=c("AB" ,"BA")),
                r299     = list("R"=c("OMPzzc" ,"UMPzzc"), "C"=c(-1 ,1),     "E"="v299",   "T"=c("AB" ,"AB")),
                r300     = list("R"=c("ATPzzc" ,"ADPzzc" ,"UTPzzc" ,"Glnzzc" ,"Gluzzc" ,"CTPzzc"), "C"=c(-1 ,1 ,-1 ,-1 ,1 ,1),     "E"="v300",   "T"=c("ABCDE" ,"ABCDE" ,"HI" ,"FG" ,"F" ,"GHI")),
                r301     = list("R"=c("ATPzzc" ,"ADPzzc" ,"AMPzzc"), "C"=c(-1 ,2 ,-1),     "E"="v301",   "T"=c("FGHIJ" ,"ABCDE" ,"ABCDE")),
                r303     = list("R"=c("ATPzzc" ,"ADPzzc" ,"GMPzzc" ,"GDPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v303",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJ" ,"FGHIJ")),
                r304     = list("R"=c("ATPzzc" ,"ADPzzc" ,"UDPzzc" ,"UMPzzc"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v304",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"FG")),
                r305     = list("R"=c("CTPzzc" ,"CDPzzc"), "C"=c(-1 ,1),     "E"="v305",   "T"=c("ABC" ,"ABC")),
                r306     = list("R"=c("ADPzzc" ,"dADPzzc"), "C"=c(-1 ,1),     "E"="v306",   "T"=c("ABCDE" ,"ABCDE")),
                r307     = list("R"=c("CDPzzc" ,"dCDPzzc"), "C"=c(-1 ,1),     "E"="v307",   "T"=c("ABC" ,"ABC")),
                r308     = list("R"=c("GDPzzc" ,"dGDPzzc"), "C"=c(-1 ,1),     "E"="v308",   "T"=c("ABCDE" ,"ABCDE")),
                r309     = list("R"=c("UDPzzc" ,"dUDPzzc"), "C"=c(-1 ,1),     "E"="v309",   "T"=c("AB" ,"AB")),
                r310     = list("R"=c("dUTPzzc" ,"dUMPzzc"), "C"=c(-1 ,1),     "E"="v310",   "T"=c("AB" ,"AB")),
                r311     = list("R"=c("DHFzzc" ,"MxxTHFzzc" ,"dUMPzzc" ,"dTMPzzc"), "C"=c(1 ,-1 ,-1 ,1),     "E"="v311",   "T"=c("CHDIEFG" ,"CDEFGHI" ,"AB" ,"AB")),
                r312     = list("R"=c("ATPzzc" ,"ADPzzc" ,"dTMPzzc" ,"dTDPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v312",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"FG")),
                r313     = list("R"=c("ATPzzc" ,"ADPzzc" ,"GDPzzc" ,"GTPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v313",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJ" ,"FGHIJ")),
                r314     = list("R"=c("ADPzzh" ,"ATPzzh" ,"GTPzzh" ,"GDPzzh"), "C"=c(1 ,-1 ,1 ,-1),     "E"="v314",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJ" ,"FGHIJ")),
                r315     = list("R"=c("ATPzzc" ,"ADPzzc" ,"UTPzzc" ,"UDPzzc"), "C"=c(-1 ,1 ,1 ,-1),     "E"="v315",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"FG")),
                r316     = list("R"=c("ATPzzc" ,"ADPzzc" ,"dADPzzc" ,"dATPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v316",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJ" ,"FGHIJ")),
                r317     = list("R"=c("ATPzzc" ,"ADPzzc" ,"dCDPzzc" ,"dCTPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v317",   "T"=c("ABCDE" ,"ABCDE" ,"FGH" ,"FGH")),
                r318     = list("R"=c("ATPzzc" ,"ADPzzc" ,"dGDPzzc" ,"dGTPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v318",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJ" ,"FGHIJ")),
                r319     = list("R"=c("ATPzzc" ,"ADPzzc" ,"dTDPzzc" ,"dTTPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v319",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"FG")),
                r320     = list("R"=c("ATPzzc" ,"ADPzzc" ,"dUDPzzc" ,"dUTPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v320",   "T"=c("ABCDE" ,"ABCDE" ,"FG" ,"FG")),
                r323     = list("R"=c("ADPzzh" ,"ATPzzh" ,"ATPzzc" ,"ADPzzc"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v323",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJ" ,"FGHIJ")),
                r326     = list("R"=c("ADPzzh" ,"ADPzzc"), "C"=c(-1 ,1),     "E"="v326",   "T"=c("ABCDE" ,"ABCDE")),
                r327     = list("R"=c("ATPzzh" ,"ATPzzc"), "C"=c(1 ,-1),     "E"="v327",   "T"=c("ABCDE" ,"ABCDE")),
                r328     = list("R"=c("NADzzc" ,"NADzzh" ,"AMPzzh" ,"AMPzzc"), "C"=c(1 ,-1 ,1 ,-1),     "E"="v328",   "T"=c("FGHIJKL" ,"FGHIJKL" ,"ABCDE" ,"ABCDE")),
                r329     = list("R"=c("ADPzzh" ,"ADPzzc" ,"NADzzc" ,"NADzzh"), "C"=c(-1 ,1 ,-1 ,1),     "E"="v329",   "T"=c("ABCDE" ,"ABCDE" ,"FGHIJKL" ,"FGHIJKL")),
                r343     = list("R"=c("Gluzzh" ,"Gluzzc"), "C"=c(-1 ,1),     "E"="v343",   "T"=c("A" ,"A")),
                r357     = list("R"=c("Alazzm" ,"Alazzc"), "C"=c(-1 ,1),     "E"="v357",   "T"=c("A" ,"A")),
                r358     = list("R"=c("Argzzm" ,"Argzzc"), "C"=c(1 ,-1),     "E"="v358",   "T"=c("ABCD" ,"ABCD")),
                r359     = list("R"=c("Gluzzm" ,"Gluzzc"), "C"=c(1 ,-1),     "E"="v359",   "T"=c("A" ,"A")),
                r364     = list("R"=c("Prozzc" ,"Prozzm"), "C"=c(-1 ,1),     "E"="v364",   "T"=c("A" ,"A")),
                r367     = list("R"=c("P5Czzc" ,"P5Czzm"), "C"=c(1 ,-1),     "E"="v367",   "T"=c("A" ,"A")),
                r370     = list("R"=c("HxxCyszzh" ,"HxxCyszzc"), "C"=c(-1 ,1),     "E"="v370",   "T"=c("A" ,"A")),
                r371     = list("R"=c("OROzzm" ,"OROzzc"), "C"=c(-1 ,1),     "E"="v371",   "T"=c("AB" ,"AB")),
                r372     = list("R"=c("DHOzzh" ,"DHOzzc"), "C"=c(-1 ,1),     "E"="v372",   "T"=c("AB" ,"AB")),
                r373     = list("R"=c("DHOzzm" ,"DHOzzc"), "C"=c(1 ,-1),     "E"="v373",   "T"=c("AB" ,"AB")),
                r374     = list("R"=c("Glyzzm" ,"Serzzm" ,"Serzzc" ,"Glyzzc"), "C"=c(1 ,-1 ,1 ,-1),     "E"="v374",   "T"=c("A" ,"B" ,"B" ,"A")),
                r384     = list("R"=c("NH4zzh" ,"NH4zzc"), "C"=c(1 ,-1),     "E"="v384",   "T"=c("A" ,"A")),
                r385     = list("R"=c("NH4zzm" ,"NH4zzc"), "C"=c(-1 ,1),     "E"="v385",   "T"=c("A" ,"A")),
                r398     = list("R"=c("Alazzp" ,"Alazzc"), "C"=c(-1 ,1),     "E"="v398",   "T"=c("A" ,"A")),
                r403     = list("R"=c("Gluzzp" ,"Gluzzc"), "C"=c(1 ,-1),     "E"="v403",   "T"=c("A" ,"A")),
                r419     = list("R"=c("Argzzh" ,"Argzzc"), "C"=c(-1 ,1),     "E"="v419",   "T"=c("ABCD" ,"ABCD")),
                r422     = list("R"=c("Gluzzh" ,"Gluzzc"), "C"=c(1 ,-1),     "E"="v422",   "T"=c("A" ,"A")),
                r423     = list("R"=c("Hiszzh" ,"Hiszzc"), "C"=c(-1 ,1),     "E"="v423",   "T"=c("ABC" ,"ABC")),
                r425     = list("R"=c("Ilezzh" ,"Ilezzc"), "C"=c(-1 ,1),     "E"="v425",   "T"=c("A" ,"A")),
                r427     = list("R"=c("Leuzzh" ,"Leuzzc"), "C"=c(-1 ,1),     "E"="v427",   "T"=c("A" ,"A")),
                r429     = list("R"=c("Lyszzh" ,"Lyszzc"), "C"=c(-1 ,1),     "E"="v429",   "T"=c("AB" ,"AB")),
                r432     = list("R"=c("Phezzh" ,"Phezzc"), "C"=c(-1 ,1),     "E"="v432",   "T"=c("A" ,"A")),
                r436     = list("R"=c("Serzzm" ,"Serzzc"), "C"=c(1 ,-1),     "E"="v436",   "T"=c("A" ,"A")),
                r437     = list("R"=c("Thrzzc" ,"Thrzzh"), "C"=c(1 ,-1),     "E"="v437",   "T"=c("A" ,"A")),
                r439     = list("R"=c("Trpzzh" ,"Trpzzc"), "C"=c(-1 ,1),     "E"="v439",   "T"=c("AB" ,"AB")),
                r441     = list("R"=c("Tyrzzh" ,"Tyrzzc"), "C"=c(-1 ,1),     "E"="v441",   "T"=c("A" ,"A")),
                r443     = list("R"=c("Valzzh" ,"Valzzc"), "C"=c(-1 ,1),     "E"="v443",   "T"=c("A" ,"A")),
                r449     = list("R"=c("ATPzzc" ,"ADPzzc"), "C"=c(-3 ,3),     "E"="v449",   "T"=c("ABCDE" ,"ABCDE")),
                r451     = list("R"=c("NH4zzx", "NH4zzc"), "C"=c(-1,1),     "E"="v451",   "T"=c("A","A")),
                r452     = list("R"=c("ATPzzc" ,"ADPzzc"), "C"=c(-3 ,3),     "E"="v452",   "T"=c("ABCDE" ,"ABCDE")),
                r543     = list("R"=c("ATPzzc" ,"ADPzzc"), "C"=c(-1 ,1),     "E"="v543",   "T"=c("ABCDE" ,"ABCDE")),
                r545     = list("R"=c("ADPzzm" ,"ATPzzm"), "C"=c(1 ,-1),     "E"="v545",   "T"=c("ABCDE" ,"ABCDE")),
                r549x5     = list("R"=c("ATPzzc"), "C"=c(-0.58559),     "E"="v549x5",   "T"=c("ABCDE")),
                r549x6     = list("R"=c("UTPzzc"), "C"=c(-0.61307),     "E"="v549x6",   "T"=c("AB")),
                r549x17     = list("R"=c("Ornzzh"), "C"=c(-1.5865),     "E"="v549x17",   "T"=c("AB")),
                r549x18     = list("R"=c("Glnzzc"), "C"=c(-51.2643),     "E"="v549x18",   "T"=c("AB")),
                r549x19     = list("R"=c("Aspzzc"), "C"=c(-101.038),     "E"="v549x19",   "T"=c("A")),
                r549x20     = list("R"=c("Gluzzc"), "C"=c(-116.394),     "E"="v549x20",   "T"=c("A")),
                r549x21     = list("R"=c("Asnzzc"), "C"=c(-70.5431),     "E"="v549x21",   "T"=c("AB")),
                r549x22     = list("R"=c("Serzzc"), "C"=c(-147.4985),     "E"="v549x22",   "T"=c("A")),
                r549x23     = list("R"=c("Cyszzc"), "C"=c(-31.1176),     "E"="v549x23",   "T"=c("A")),
                r549x24     = list("R"=c("GABAzzc"), "C"=c(-0.09401),     "E"="v549x24",   "T"=c("A")),
                r549x25     = list("R"=c("Thrzzc"), "C"=c(-136.6831),     "E"="v549x25",   "T"=c("A")),
                r549x26     = list("R"=c("Glyzzc"), "C"=c(-273.4827),     "E"="v549x26",   "T"=c("A")),
                r549x27     = list("R"=c("Metzzc"), "C"=c(-33.615),     "E"="v549x27",   "T"=c("A")),
                r549x28     = list("R"=c("ureazzm"), "C"=c(-2.8404),     "E"="v549x28",   "T"=c("AB")),
                r549x29     = list("R"=c("Prozzc"), "C"=c(-108.0369),     "E"="v549x29",   "T"=c("A")),
                r549x30     = list("R"=c("CTPzzc"), "C"=c(-0.49638),     "E"="v549x30",   "T"=c("ABC")),
                r549x31     = list("R"=c("GTPzzc"), "C"=c(-0.40478),     "E"="v549x31",   "T"=c("ABCDE")),
                r549x32     = list("R"=c("dATPzzc"), "C"=c(-0.68091),     "E"="v549x32",   "T"=c("ABCDE")),
                r549x33     = list("R"=c("dCTPzzc"), "C"=c(-0.40366),     "E"="v549x33",   "T"=c("ABC")),
                r549x34     = list("R"=c("dGTPzzc"), "C"=c(-0.37141),     "E"="v549x34",   "T"=c("ABCDE")),
                r549x35     = list("R"=c("dTTPzzc"), "C"=c(-0.69281),     "E"="v549x35",   "T"=c("AB")),
                r549x36     = list("R"=c("Alazzc"), "C"=c(-230.671),     "E"="v549x36",   "T"=c("A")),
                r549x37     = list("R"=c("Argzzc"), "C"=c(-74.6825),     "E"="v549x37",   "T"=c("ABCD")),
                r549x38     = list("R"=c("Lyszzc"), "C"=c(-103.5935),     "E"="v549x38",   "T"=c("AB")),
                r549x39     = list("R"=c("Hiszzc"), "C"=c(-33.2443),     "E"="v549x39",   "T"=c("ABC")),
                r549x40     = list("R"=c("Ilezzc"), "C"=c(-89.8269),     "E"="v549x40",   "T"=c("A")),
                r549x41     = list("R"=c("Leuzzc"), "C"=c(-158.7712),     "E"="v549x41",   "T"=c("A")),
                r549x42     = list("R"=c("Phezzc"), "C"=c(-66.6617),     "E"="v549x42",   "T"=c("A")),
                r549x43     = list("R"=c("Trpzzc"), "C"=c(-16.871),     "E"="v549x43",   "T"=c("AB")),
                r549x44     = list("R"=c("Tyrzzc"), "C"=c(-47.0292),     "E"="v549x44",   "T"=c("A")),
                r549x45     = list("R"=c("Valzzc"), "C"=c(-146.2688),     "E"="v549x45",   "T"=c("A"))
    )
    
    # equations for determined fluxes
    .GlobalEnv$eq_det <- c("v178b = v178")
    
    # construct the model
    .GlobalEnv$net <- net2mat(rxn, add_eq=eq_det)
    
    ####################################
    cat("\n   ... simulate labeling dynamics ...\n\n")
    ####################################
    
    # labeling dynamics of label input (here we simulate a switch from unlabeled to fully-labeled nutrient Sout)
    .GlobalEnv$anFun <- list("NH4zzx_1-M0" = "0.0+0.0*t", "NH4zzx_1-M1" = "1.0+0.0*t")
    
    # fluxes and metabolite concentrations
    .GlobalEnv$fluxes <- c("v5"=17.8571,
                "v7"=38.7663,
                "v18"=15.8975,
                "v21"=0.47596,
                "v22"=0.23589,
                "v23"=0.23798,
                "v24"=0.0020877,
                "v41"=0.622,
                "v42"=0.034163,
                "v44"=0.0051551,
                "v47"=0.00046189,
                "v49"=0.29111,
                "v50"=0.29111,
                "v52"=9.7314,
                "v57"=9.1055,
                "v58"=9.1055,
                "v59"=0.70108,
                "v62"=38.7663,
                "v67"=0.41566,
                "v68"=0.41566,
                "v69"=0.41566,
                "v70"=0.41566,
                "v71"=0.41566,
                "v75"=7.7504e-13,
                "v76"=7.7504e-13,
                "v77"=7.7504e-13,
                "v80"=8.5858,
                "v81"=9.5515,
                "v84"=23.962,
                "v90"=0.26738,
                "v91"=0.26738,
                "v92"=0.26738,
                "v93"=0.26738,
                "v100"=37.4649,
                "v101"=37.4649,
                "v106"=0.10475,
                "v107"=0.10445,
                "v110"=0.043492,
                "v111"=0.96505,
                "v113"=0.14052,
                "v114"=1.1107,
                "v124"=0.28165,
                "v129"=0.054556,
                "v130"=0.07945,
                "v131"=0.77789,
                "v132"=0.62398,
                "v133"=0.62398,
                "v137"=0.029864,
                "v138"=0.029864,
                "v139"=0.00055425,
                "v140"=0.001719,
                "v148"=0.051787,
                "v149"=0.051787,
                "v150"=0.051787,
                "v155"=6.0126e-28,
                "v156"=0.18454,
                "v157"=0.063288,
                "v158"=0.063288,
                "v159"=0.063288,
                "v160"=0.063288,
                "v161"=0.063784,
                "v162"=0.062019,
                "v163"=0.062019,
                "v164"=0.062019,
                "v166"=0.056435,
                "v167"=0.13727,
                "v168"=1.1107,
                "v172"=0.024894,
                "v173"=0.026892,
                "v175"=0.024894,
                "v176"=0.026892,
                "v178"=39.8566,
                "v179"=0.69881,
                "v182"=41.589,
                "v184"=0.084158,
                "v185"=0.013468,
                "v186"=7.5209e-05,
                "v187"=3.2384e-29,
                "v188"=3.2383e-29,
                "v189"=1.278e-33,
                "v190"=0.098788,
                "v191"=1.8395,
                "v193"=0.72611,
                "v194"=0.72611,
                "v199"=0.0016342,
                "v200"=0.026596,
                "v201"=0.026596,
                "v202"=0.026596,
                "v203"=0.026596,
                "v204"=0.026596,
                "v205"=0.026596,
                "v206"=0.026596,
                "v207"=0.026596,
                "v208"=0.026596,
                "v209"=0.026596,
                "v210"=0.071862,
                "v211"=0.3159,
                "v212"=0.071862,
                "v214"=0.071862,
                "v216"=0.071862,
                "v217"=0.24403,
                "v219"=0.24403,
                "v221"=0.12702,
                "v224"=0.12702,
                "v226"=0.12702,
                "v227"=1.0171,
                "v228"=1.0171,
                "v229"=0.082876,
                "v231"=0.082876,
                "v232"=0.082876,
                "v233"=0.082876,
                "v234"=0.082876,
                "v235"=0.93421,
                "v237"=0.93421,
                "v238"=0.026892,
                "v239"=0.026892,
                "v240"=0.026892,
                "v241"=0.026892,
                "v242"=0.027446,
                "v247"=0.090954,
                "v248"=0.05333,
                "v249"=0.0022724,
                "v251"=0.0022724,
                "v252"=1.2133e-11,
                "v254"=0.084158,
                "v255"=1.2133e-11,
                "v257"=0.084158,
                "v258"=1.2133e-11,
                "v260"=0.08643,
                "v261"=0.35558,
                "v262"=0.26915,
                "v267"=0.040389,
                "v268"=0.040389,
                "v269"=0.040389,
                "v270"=0.097041,
                "v271"=0.90732,
                "v272"=0.013497,
                "v273"=0.013497,
                "v274"=0.013497,
                "v275"=0.013497,
                "v276"=0.013497,
                "v277"=0.013497,
                "v278"=0.037624,
                "v279"=0.11702,
                "v280"=0.0016342,
                "v281"=0.0016342,
                "v282"=0.0016342,
                "v283"=0.0016342,
                "v284"=0.0016342,
                "v285"=0.0016342,
                "v286"=0.0016342,
                "v287"=0.0016342,
                "v288"=0.02823,
                "v289"=0.02823,
                "v290"=0.00062096,
                "v291"=0.00062096,
                "v292"=0.00062096,
                "v293"=0.02823,
                "v294"=0.02823,
                "v295"=0.0017648,
                "v296"=0.0017648,
                "v297"=0.0017648,
                "v298"=0.0017648,
                "v299"=0.0017648,
                "v300"=0.00072004,
                "v301"=0.054801,
                "v303"=0.00062096,
                "v304"=0.0017648,
                "v305"=0.00032293,
                "v306"=0.00054473,
                "v307"=0.00032293,
                "v308"=0.00029713,
                "v309"=0.00055425,
                "v310"=0.00055425,
                "v311"=0.00055425,
                "v312"=0.00055425,
                "v313"=0.00032383,
                "v314"=0.02823,
                "v315"=0.62321,
                "v316"=0.00054473,
                "v317"=0.00032293,
                "v318"=0.00029713,
                "v319"=0.00055425,
                "v320"=0.00055425,
                "v323"=6.3399,
                "v326"=0.0016342,
                "v327"=1.0651e-28,
                "v328"=0.0016342,
                "v329"=0.0016342,
                "v343"=1.2512,
                "v357"=5.6884e-28,
                "v358"=0.0022724,
                "v359"=0.17951,
                "v364"=0.26915,
                "v367"=0.35558,
                "v370"=0.026892,
                "v371"=0.0017648,
                "v372"=0.0017648,
                "v373"=0.0017648,
                "v374"=0.53476,
                "v384"=1.7408,
                "v385"=0.36501,
                "v398"=0.18454,
                "v403"=0.18454,
                "v419"=0.062019,
                "v422"=1.3145,
                "v423"=0.026596,
                "v425"=0.071862,
                "v427"=0.12702,
                "v429"=0.082876,
                "v432"=0.05333,
                "v436"=0.26738,
                "v437"=0.83545,
                "v439"=0.013497,
                "v441"=0.037624,
                "v443"=0.11702,
                "v449"=0.010197,
                "v451"=2.0757,
                "v452"=0.051787,
                "v543"=1.6597,
                "v545"=23.8778,
                "v549x5"=0.00080001,
                "v549x6"=0.00080001,
                "v549x17"=0.00080001,
                "v549x18"=0.00080001,
                "v549x19"=0.00080001,
                "v549x20"=0.00080001,
                "v549x21"=0.00080001,
                "v549x22"=0.00080001,
                "v549x23"=0.00080001,
                "v549x24"=0.00080001,
                "v549x25"=0.00080001,
                "v549x26"=0.00080001,
                "v549x27"=0.00080001,
                "v549x28"=0.00080001,
                "v549x29"=0.00080001,
                "v549x30"=0.00080001,
                "v549x31"=0.00080001,
                "v549x32"=0.00080001,
                "v549x33"=0.00080001,
                "v549x34"=0.00080001,
                "v549x35"=0.00080001,
                "v549x36"=0.00080001,
                "v549x37"=0.00080001,
                "v549x38"=0.00080001,
                "v549x39"=0.00080001,
                "v549x40"=0.00080001,
                "v549x41"=0.00080001,
                "v549x42"=0.00080001,
                "v549x43"=0.00080001,
                "v549x44"=0.00080001,
                "v549x45"=0.00080001
    )
    
    .GlobalEnv$meta_conc <- c("NADPzzh"=1.100000e-02,
                   "NADPHzzh"=1.100000e-02,
                   "ADPzzh"=1.559359e-01,
                   "ATPzzh"=1.568179e-01,
                   "ADPGzzh"=2.200000e-02,
                   "ATPzzc"=2.345588e-02,
                   "ADPzzc"=2.369978e-02,
                   "UTPzzc"=1.100000e-01,
                   "UDPGzzc"=2.200000e+00,
                   "UDPzzc"=2.200000e-01,
                   "NADzzc"=1.100000e-02,
                   "NADHzzc"=1.100000e-02,
                   "NADPzzc"=1.100000e-02,
                   "NADPHzzc"=1.100000e-02,
                   "NADzzh"=1.100000e-02,
                   "NADHzzh"=1.100000e-02,
                   "ThPPzzm"=1.100000e-02,
                   "HxxEthxxThPPzzm"=1.100000e-02,
                   "LPAzzm"=1.100000e-02,
                   "AxxDHLzzm"=1.100000e-02,
                   "CoAzzm"=1.100000e-02,
                   "AxxCoAzzm"=1.100000e-02,
                   "DHLzzm"=1.100000e-02,
                   "NADzzm"=1.100000e-02,
                   "NADHzzm"=1.100000e-02,
                   "SxxDHLzzm"=1.100000e-02,
                   "SxxCoAzzm"=1.100000e-02,
                   "ADPzzm"=6.236434e-02,
                   "ATPzzm"=6.172623e-02,
                   "Gluzzp"=8.168819e-02,
                   "Glyzzm"=5.016197e+00,
                   "LPLzzm"=1.100000e-02,
                   "amDHPzzm"=1.100000e-01,
                   "THFzzm"=1.100000e-02,
                   "DHPzzm"=1.100000e-02,
                   "MxxTHFzzm"=1.100000e-02,
                   "NH4zzm"=2.116383e+00,
                   "Serzzm"=1.483960e+01,
                   "AMPzzh"=1.038596e-01,
                   "NADPHzzm"=1.100000e-02,
                   "NADPzzm"=1.100000e-02,
                   "AxxCoAzzc"=1.100000e-02,
                   "CoAzzc"=1.100000e-02,
                   "AxxCoAzzh"=1.100000e-02,
                   "CoAzzh"=1.100000e-02,
                   "MxxCoAzzh"=1.100000e-02,
                   "THFzzh"=1.100000e-02,
                   "FxxTHFzzh"=1.100000e-02,
                   "DHFzzc"=1.100000e-02,
                   "THFzzc"=1.100000e-02,
                   "AMPzzc"=1.100000e-02,
                   "NH4zzh"=1.066581e+01,
                   "APSzzh"=1.100000e-01,
                   "GSHzzh"=1.100000e-02,
                   "GSSGzzh"=1.100000e-02,
                   "Gluzzh"=3.744905e+01,
                   "Gluzzm"=8.046879e-02,
                   "Alazzm"=1.100000e-02,
                   "Alazzp"=3.300000e+00,
                   "AxxGluzzh"=4.620000e-01,
                   "AxxGluPzzh"=1.100000e-01,
                   "AxxGluxxSeAzzh"=1.100000e-01,
                   "AxxOrnzzh"=1.540000e-01,
                   "Ornzzh"=1.486623e+00,
                   "Glnzzh"=9.876600e+01,
                   "CBPzzh"=1.100000e-01,
                   "CTLzzh"=1.760000e+01,
                   "Aspzzh"=1.223765e+01,
                   "ArgxxSCAzzh"=1.100000e+00,
                   "Argzzh"=7.561474e+00,
                   "Glnzzc"=2.339996e-01,
                   "Aspzzc"=1.512351e+00,
                   "Gluzzc"=8.887904e-01,
                   "Asnzzc"=1.650000e+01,
                   "NH4zzc"=1.471781e+01,
                   "Serzzc"=1.483960e+01,
                   "AxxSerzzc"=1.100000e-02,
                   "Serzzh"=1.120792e+00,
                   "AxxSerzzh"=1.100000e-02,
                   "Cyszzc"=1.057562e-01,
                   "Cyszzh"=1.142438e-01,
                   "GABAzzc"=1.100000e-01,
                   "GABAzzm"=1.100000e-02,
                   "Thrzzc"=3.691247e+00,
                   "Glyzzc"=7.068475e+00,
                   "Glyzzh"=1.532882e-02,
                   "PRxxATPzzh"=1.100000e-01,
                   "PRxxAMPzzh"=1.100000e-01,
                   "PxxAICARxxPzzh"=1.100000e-01,
                   "PuxxAICARxxPzzh"=1.100000e-01,
                   "EIGPzzh"=1.100000e-01,
                   "AICARzzh"=1.100000e-01,
                   "IAxxPzzh"=1.100000e-01,
                   "HisolxxPzzh"=1.100000e-01,
                   "Hisolzzh"=1.100000e-01,
                   "Hisalzzh"=1.100000e-01,
                   "Hiszzh"=2.200000e-01,
                   "Thrzzh"=4.008753e+00,
                   "ThPPzzh"=1.100000e-02,
                   "HxxEthxxThPPzzh"=1.100000e-02,
                   "Ilezzh"=2.750000e-01,
                   "Leuzzh"=2.750000e-01,
                   "AspPzzh"=1.100000e-01,
                   "AspxxSeAzzh"=1.100000e-01,
                   "DHDzzh"=1.100000e-01,
                   "THDzzh"=1.100000e-01,
                   "DAPzzh"=1.100000e-01,
                   "mDAPzzh"=1.100000e-01,
                   "Lyszzh"=1.650000e-01,
                   "HxxSerzzh"=1.100000e-01,
                   "PHxxSerzzh"=1.100000e-01,
                   "CTHzzh"=1.100000e-01,
                   "HxxCyszzh"=5.500000e-02,
                   "MxxTHFzzc"=1.100000e-02,
                   "5MxxTHFzzc"=1.100000e-02,
                   "HxxCyszzc"=5.500000e-02,
                   "Metzzc"=2.750000e-01,
                   "AGNzzh"=1.100000e-01,
                   "Phezzh"=3.850000e-01,
                   "Argzzm"=2.770515e-01,
                   "Ornzzm"=5.337729e-02,
                   "ureazzm"=2.750000e+00,
                   "GluxxSeAzzm"=1.100000e+00,
                   "GluPzzc"=1.100000e-02,
                   "GluPzzm"=1.100000e+00,
                   "GluxxSeAzzc"=1.100000e-02,
                   "P5Czzc"=5.500000e-02,
                   "P5Czzm"=5.500000e-02,
                   "Prozzc"=1.126966e+00,
                   "Prozzm"=8.530339e-01,
                   "PSerzzh"=1.100000e-01,
                   "ANTzzh"=1.100000e-01,
                   "PRxxANTzzh"=1.100000e-01,
                   "CPDxxRu5Pzzh"=1.100000e-01,
                   "IndxxGPzzh"=1.100000e-01,
                   "Indzzh"=1.100000e-01,
                   "Trpzzh"=1.650000e-01,
                   "Tyrzzh"=2.200000e-01,
                   "Valzzh"=5.500000e-01,
                   "PRAzzh"=1.100000e-01,
                   "GARzzh"=1.100000e-01,
                   "FGARzzh"=1.100000e-01,
                   "FGAMzzh"=1.100000e-01,
                   "AIRzzh"=1.100000e-01,
                   "CAIRzzh"=1.100000e-01,
                   "SAICARzzh"=1.100000e-01,
                   "FAICARzzh"=1.100000e-01,
                   "IMPzzh"=1.076325e-01,
                   "IMPzzc"=1.100000e-02,
                   "XMPzzc"=1.100000e-01,
                   "GMPzzc"=1.100000e-01,
                   "GTPzzh"=1.740040e-01,
                   "DCxxAMPzzh"=1.100000e-01,
                   "GDPzzh"=1.076325e-01,
                   "CBxxAspzzh"=1.100000e-01,
                   "DHOzzh"=3.666667e-02,
                   "DHOzzm"=3.666667e-02,
                   "OROzzm"=5.500000e-02,
                   "OROzzc"=5.500000e-02,
                   "OMPzzc"=1.100000e-01,
                   "UMPzzc"=4.070000e+00,
                   "CTPzzc"=1.100000e-01,
                   "GDPzzc"=1.100000e-02,
                   "CDPzzc"=1.100000e-01,
                   "dADPzzc"=1.100000e-01,
                   "dCDPzzc"=1.100000e-01,
                   "dGDPzzc"=1.100000e-01,
                   "dUDPzzc"=1.100000e-01,
                   "dUTPzzc"=1.100000e-01,
                   "dUMPzzc"=1.100000e-01,
                   "dTMPzzc"=1.100000e-01,
                   "dTDPzzc"=1.100000e-01,
                   "GTPzzc"=1.100000e-02,
                   "dATPzzc"=1.100000e-01,
                   "dCTPzzc"=1.100000e-01,
                   "dGTPzzc"=1.100000e-01,
                   "dTTPzzc"=1.100000e-01,
                   "Alazzc"=3.300000e+00,
                   "Argzzc"=7.561474e+00,
                   "Lyszzc"=1.650000e-01,
                   "Hiszzc"=2.200000e-01,
                   "DHOzzc"=3.666667e-02,
                   "Ilezzc"=2.750000e-01,
                   "Leuzzc"=2.750000e-01,
                   "Phezzc"=3.850000e-01,
                   "Trpzzc"=1.650000e-01,
                   "Tyrzzc"=2.200000e-01,
                   "Valzzc"=5.500000e-01
    )
    .GlobalEnv$meta_conc_org = meta_conc
}

simulate_aracore_model <- function() {
    # simulation times (exponential sampling frequency from 0 to 15 min)
    #times <- round(10**(seq(0, log10(16), length.out=30))-1, 2)
    .GlobalEnv$times <- seq(0, 320, by=4);
    #times <- seq(0, 6);
    
    # run simulations
    .GlobalEnv$res <- simulate(net       = net,
                    kp        = fluxes,
                    anFun     = anFun,
                    p         = 0.0,
                    trf       = xch2fb,
                    meta_conc = meta_conc,
                    method    = "R",
                    unloadDll = TRUE,
                    times     = times)
    .GlobalEnv$res_org = res;
}

relevant_metabolites = c("AGNzzh_1" ,"Phezzh_1", "Gluzzh_1" ,"Ilezzh_1", "HxxCyszzc_1" ,"Metzzc_1","Aspzzh_1", "Thrzzc_1" ,"Glyzzc_1", "Serzzc_1");
subset_times = c(1,17,33,49,65,81)

run_fit_input <- function() {
    ####################################
    cat("\n   ... fit labeling dynamics of the relevant metabolites ...\n\n")
    ####################################
    
    # fit analytical functions to the labeling dynamics of all metabolites
    enr_input <- res$res_dyn$enrichments[subset_times,relevant_metabolites]
    
    # # fit analytical functions to the labeling dynamics of all metabolites
    .GlobalEnv$enr_in <- fit_label_input(enr_input, t=times[subset_times], file="res_fit_enr", mc.cores=numCores)
}

run_one_estimation <- function() {
    
    ####################################
    cat("\n   ... calculate fluxes for all minimal subsystems ...\n\n")
    ####################################
    
    # definition of all minimal subsystems to analyze
    # here, fluxes and pools are estimated for all minimal subsystems of the network provided as example
    subsystems <- list(
        s248 = list(name = "s248",
                    rxn_subnet = list(r248     = list("R"=c("AGNzzh" ,"Phezzh"), "C"=c(-1 ,1),     "E"="v248",   "T"=c("A" ,"A")),
                                      rout = list("R"=c("Phezzh"),        "C"=c(-1),    "E"="v248", "T"=c("A"))),
                    meta_conc_subnet = c("AGNzzh"=1, "Phezzh"=1.0),
                    kp_subnet = c("v248"=0.5),
                    te_subnet = c("v248", "Phezzh"),
                    te_upc_subnet = c("v248"=100, "Phezzh"=100),
                    te_loc_subnet = c("v248"=1e-3, "Phezzh"=0.01),
                    sd_meas = list(iso=0.02, conc=c("Phezzh"=0.05)),
                    times = times[subset_times],
                    enr_in = enr_in,
                    anFun = NULL,
                    niter = NULL,
                    mc.cores = NULL,
                    data_meas_subnet = list(conc=c("Phezzh"=0.385), iso=cbind(times[subset_times], "Phezzh_1"=res$res_dyn$enrichments[subset_times, "Phezzh_1"])))
        ,
        s216 = list(name = "s216",
                    rxn_subnet = list(r216     = list("R"=c("Gluzzh" ,"Ilezzh"), "C"=c(-1 ,1),     "E"="v216",   "T"=c("A" ,"A")),
                                      rout = list("R"=c("Ilezzh"),        "C"=c(-1),    "E"="v216", "T"=c("A"))),
                    meta_conc_subnet = c("Gluzzh"=1, "Ilezzh"=1.0),
                    kp_subnet = c("v216"=0.5),
                    te_subnet = c("v216", "Ilezzh"),
                    te_upc_subnet = c("v216"=100, "Ilezzh"=100),
                    te_loc_subnet = c("v216"=1e-3, "Ilezzh"=0.01),
                    sd_meas = list(iso=0.02, conc=c("Ilezzh"=0.05)),
                    times = times[subset_times],
                    enr_in = enr_in,
                    anFun = NULL,
                    niter = NULL,
                    mc.cores = NULL,
                    data_meas_subnet = list(conc=c("Ilezzh"=0.275), iso=cbind(times[subset_times], "Ilezzh_1"=res$res_dyn$enrichments[subset_times, "Ilezzh_1"])))
        ,
        s241 = list(name = "s241",
                    rxn_subnet = list(r241     = list("R"=c("HxxCyszzc" ,"Metzzc"), "C"=c(-1 ,1),     "E"="v241",   "T"=c("A" ,"A")),
                                      rout = list("R"=c("Metzzc"),        "C"=c(-1),    "E"="v241", "T"=c("A"))),
                    meta_conc_subnet = c("HxxCyszzc"=1, "Metzzc"=1.0),
                    kp_subnet = c("v241"=0.5),
                    te_subnet = c("v241", "Metzzc"),
                    te_upc_subnet = c("v241"=100, "Metzzc"=100),
                    te_loc_subnet = c("v241"=1e-3, "Metzzc"=0.01),
                    sd_meas = list(iso=0.02, conc=c("Metzzc"=0.5)),
                    times = times[subset_times],
                    enr_in = enr_in,
                    anFun = NULL,
                    niter = NULL,
                    mc.cores = NULL,
                    data_meas_subnet = list(conc=c("Metzzc"=0.275), iso=cbind(times[subset_times], "Metzzc_1"=res$res_dyn$enrichments[subset_times, "Metzzc_1"])))
        ,
        s168 = list(name = "s168",
                    rxn_subnet = list(r168   = list("R"=c("Gluzzh" ,"Aspzzh"), "C"=c(-1 ,1),     "E"="v168",   "T"=c("A" ,"A")),
                                      rout = list("R"=c("Aspzzh"),        "C"=c(-1),    "E"="v168", "T"=c("A"))),
                    meta_conc_subnet = c("Gluzzh"=1, "Aspzzh"=1.0),
                    kp_subnet = c("v168"=0.5),
                    te_subnet = c("v168", "Aspzzh"),
                    te_upc_subnet = c("v168"=100, "Aspzzh"=100),
                    te_loc_subnet = c("v168"=1e-3, "Aspzzh"=0.01),
                    sd_meas = list(iso=0.02, conc=c("Aspzzh"=0.5)),
                    times = times[subset_times],
                    enr_in = enr_in,
                    anFun = NULL,
                    niter = NULL,
                    mc.cores = NULL,
                    data_meas_subnet = list(conc=c("Aspzzh"=12.2), iso=cbind(times[subset_times], "Aspzzh_1"=res$res_dyn$enrichments[subset_times, "Aspzzh_1"])))
        ,
        s193_242 = list(name = "s193_242",
                    rxn_subnet = list(r193   = list("R"=c("Thrzzc" ,"Glyzzc"), "C"=c(-1 ,1),     "E"="v193",   "T"=c("A" ,"A")),
                                      r242 = list("R"=c("Serzzc" ,"Glyzzc"), "C"=c(-1 ,1),     "E"="v242",   "T"=c("A" ,"A" )),
                                      rout = list("R"=c("Glyzzc"),        "C"=c(-1),    "E"="v193+v242", "T"=c("A"))),
                    meta_conc_subnet = c("Glyzzc"=1, "Thrzzc"=1.0, "Serzzc"=1.0),
                    kp_subnet = c("v193"=0.5, "v242"=0.5),
                    te_subnet = c("v193", "v242", "Glyzzc"),
                    te_upc_subnet = c("v193" = 100, "v242"=100, "Glyzzc"=100),
                    te_loc_subnet = c("v193" = 1E-4, "v242"=1E-4, "Glyzzc"=1E-3),
                    sd_meas = list(iso=0.02, conc=c("Glyzzc"=0.5)),
                    times = times[subset_times],
                    enr_in = enr_in,
                    anFun = NULL,
                    niter = NULL,
                    mc.cores = NULL,
                    data_meas_subnet = list(conc=c("Glyzzc"=7), iso=cbind(times[subset_times], "Glyzzc_1"=res$res_dyn$enrichments[subset_times, "Glyzzc_1"])))
    )
    
    # calculate fluxes for all subsystems
    res_sub <- fit_subsystems(subsystems, dirname="fit_minimal_subsystems", mc.cores=numCores)
    return(res_sub)
}

setup_env()
init_aracore_model()
simulate_aracore_model()

mc_iter = 100
meta_conc_rel_sd = 0.05
isotopomer_abs_sd = 0.02

start_time <- Sys.time()
meta_conc_rel_noise = 1
isotopomer_abs_noise = 0
timepoint_selections = list(seq(1,81,16), seq(1,81,8), seq(1,81,4), seq(1,81))
names(timepoint_selections) = c("6", "11", "21", "81")
all_enr_in = list(list(), list(), list(), list())
names(all_enr_in) = c("6", "11", "21", "81")
all_res_iso = list(list(), list(), list(), list())
names(all_res_iso) = c("6", "11", "21", "81")
all_res_kfp = list(list(), list(), list(), list())
names(all_res_kfp) = c("6", "11", "21", "81")
all_res_ratio = list(list(), list(), list(), list())
names(all_res_ratio) = c("6", "11", "21", "81")


kfp_subsystems = list(r248=c(s="AGNzzh", p="Phezzh"), r216=c(s="Gluzzh", p="Ilezzh"), 
                      r241=c(s="HxxCyszzc", p="Metzzc"), r168=c(s="Gluzzh", p="Aspzzh"))

estimate_kfp <- function() {
  kfp_res = list()
  for (kfp_sub in names(kfp_subsystems)) {
    v_part = t(enr_in[["all"]][[paste0(kfp_subsystems[[kfp_sub]][["s"]],"_1")]]$exp)
    w_part = t(enr_in[["all"]][[paste0(kfp_subsystems[[kfp_sub]][["p"]], "_1")]]$exp)
    kfp_res[[kfp_sub]] = kfp_simple(times[subset_times], 1-v_part, 1-w_part)
  }
  return (kfp_res)
}

estimate_flux_ratio <- function() {
  fluxratio_res = list()
  fluxratio_res = fluxratio_simple3(times = times[subset_times], a_t = t(enr_in[["all"]][["Thrzzc_1"]]$exp), 
                                    b_t = t(enr_in[["all"]][["Serzzc_1"]]$exp), z_t = t(enr_in[["all"]][["Glyzzc_1"]]$exp),
                                    formula_z = enr_in$all[["Glyzzc_1"]]$anFun[["Glyzzc_1-M1"]], pz = meta_conc_org[["Glyzzc"]])
  return(fluxratio_res)
}

for (mc_run in seq(mc_iter+1)) {
    meta_conc = meta_conc_org * meta_conc_rel_noise
    res$res_dyn$enrichments <- res_org$res_dyn$enrichments + isotopomer_abs_noise
    res$res_dyn$enrichments[res$res_dyn$enrichments <0] <- 0
    res$res_dyn$enrichments[res$res_dyn$enrichments >1] <- 1
    for (name in names(timepoint_selections)) {
        .GlobalEnv$subset_times = timepoint_selections[[name]]
        print(paste0("Run #",as.character(mc_run), " for  " ,name, " timepoints"))
        run_fit_input()
        all_enr_in[[name]][[as.character(mc_run)]] = enr_in
        res_iso = run_one_estimation()
        all_res_iso[[name]][[as.character(mc_run)]] = res_iso
        res_kfp = estimate_kfp()
        all_res_kfp[[name]][[as.character(mc_run)]] = res_kfp
        res_ratio = estimate_flux_ratio()
        all_res_ratio[[name]][[as.character(mc_run)]] = res_ratio
    }
    # set noise here, so the first run is without noise
    meta_conc_rel_noise = rnorm(length(as.numeric(meta_conc)), mean=1.0, sd=meta_conc_rel_sd)
    isotopomer_abs_noise = rnorm(length(as.numeric(res$res_dyn$enrichments)), mean=0, sd=isotopomer_abs_sd)
}
end_time <- Sys.time()
end_time - start_time

all_res_kfp = list(list(), list(), list(), list())
names(all_res_kfp) = c("6", "11", "21", "81")

for (mc_run in seq(mc_iter+1)) {
    for (name in names(all_res_kfp)) {
        .GlobalEnv$subset_times = timepoint_selections[[name]]
        print(paste0("Run #",as.character(mc_run), " for  " ,name, " timepoints"))
        enr_in = all_enr_in[[name]][[as.character(mc_run)]]
        res_kfp = estimate_kfp()
        all_res_kfp[[name]][[as.character(mc_run)]] = res_kfp
    }
}


all_res_ratio = list(list(), list(), list(), list())
names(all_res_ratio) = c("6", "11", "21", "81")

for (mc_run in seq(mc_iter+1)) {
    for (name in names(all_res_ratio)) {
        .GlobalEnv$subset_times = timepoint_selections[[name]]
        print(paste0("Run #",as.character(mc_run), " for  " ,name, " timepoints"))
        enr_in = all_enr_in[[name]][[as.character(mc_run)]]
        res_ratio = estimate_flux_ratio()
        all_res_ratio[[name]][[as.character(mc_run)]] = res_ratio
    }
}

fluxratio_estimate_data = list()
fluxratio_r_position = 0
for (fluxratio_r in c("r193_", "r242_")) {
    fluxratio_r_position = 1 + fluxratio_r_position
    for (tp in c("6", "11", "21", "81")) {
        fluxratio_estimate_data[[paste0(fluxratio_r,tp)]] = c()
        for (i in seq(mc_iter+1)) {
            fluxratio_estimate_data[[paste0(fluxratio_r,tp)]][i] = all_res_ratio[[tp]][[as.character(i)]]$m$getPars()[fluxratio_r_position]
        }
    }
}



# go back to the initial working directory
setwd(wd)


####################################
cat("\nDone.\n")
####################################

flux_names = c("v168", "v193", "v216", "v241", "v242", "v248")
subsystem_names = c("s168", "s193_242", "s216", "s241", "s193_242", "s248")
estimate_data = list()
p_value_data = list()
extract_data <- function () {
    for (fn_i in seq(length(flux_names))) { 
        for (tp in c("6", "11", "21", "81")) {
            fn = flux_names[fn_i]
            sn = subsystem_names[fn_i]
            fn_tp = paste0(fn,"_",tp)
            estimate_data[[fn_tp]] = c()
            for (i in seq(mc_iter+1)) {
                estimate_data[[fn_tp]][i] = all_res_iso_0_100[[tp]][[as.character(i)]][["res"]][[sn]][["result"]][["par"]][[fn]]
                p_value_data[[fn_tp]][i] = all_res_iso_0_100[[tp]][[as.character(i)]][["res"]][[sn]][["chi2"]][["p-value, i.e. P(X^2<=value)"]]
            }
        }
    }
}
extract_data

kfp_estimate_data = list()
#extract_kfp_data <- function () {
for (kfp_sub in names(kfp_subsystems)) {
    for (tp in c("6", "11", "21", "81")) {
        kfp_sub_tp = paste0(kfp_sub,"_",tp)
        print(kfp_sub_tp)
        kfp_estimate_data[[kfp_sub_tp]] = c()
        for (i in seq(1,51)) {
            kfp_estimate_data[[kfp_sub_tp]][i] = NA
            tryCatch({ kfp_estimate_data[[kfp_sub_tp]][i] = all_res_kfp[[tp]][[as.character(i)]][[kfp_sub]][["k_w"]]*meta_conc[kfp_subsystems[[kfp_sub]][["p"]]] },
                     error=function(cond) {
                         message(cond)
                         # Choose a return value in case of error
                         return(NA)
                     },
                     warning=function(cond) {
                         message(cond)
                         # Choose a return value in case of warning
                         return(NULL)
                     })
        }
    }
}



log_estimate_data = list();
for (fn_tp in names(estimate_data)) { log_estimate_data[[fn_tp]] = log(estimate_data[[fn_tp]],10) }

log_kfp_estimate_data = list();
for (fn_tp in names(kfp_estimate_data)) { log_kfp_estimate_data[[fn_tp]] = log(kfp_estimate_data[[fn_tp]],10) }

conf_interval_data = list();
for (fn_tp in names(estimate_data)) { conf_interval_data[[fn_tp]] = quantile(estimate_data[[fn_tp]],probs=c(0.025,0.975)) }

conf_interval_diff = list();
for (fn_tp in names(conf_interval_data)) { conf_interval_diff[[fn_tp]] = conf_interval_data[[fn_tp]][2]-conf_interval_data[[fn_tp]][2] }
conf_interval_rel = list();
for (fn_tp in names(conf_interval_data)) { conf_interval_rel[[fn_tp]] = conf_interval_data[[fn_tp]][2]/conf_interval_data[[fn_tp]][1] }
conf_interval_rel_log = list();
 for (fn_tp in names(conf_interval_data)) { conf_interval_rel_log[[fn_tp]] = log(conf_interval_data[[fn_tp]][2]/conf_interval_data[[fn_tp]][1],10) }


plot_figure5 <- function() {
    pdf(file="Fig_5.pdf", width=10, height=6)
    colors = rep(fun_col(length(flux_names)), each=4)
    for (c_i in seq(length(colors))) {
        colors[c_i] = paste0(colors[c_i], as.hexmode(80+((c_i-1) %% 4)*48))
    }
    vioplot(log_estimate_data,col = colors)
    true_values = rep(log(fluxes[flux_names],10), each=4)
    points(true_values, pch=4, col="black", cex=3)
    dev.off()
    
    illus_plots = list()
    colors = rep("#AAAAAA",each=4)
    for (p_i in seq(4)) {
        illus_plots[[p_i]] = c(1,0.7,0.5,0.5,0.5,0.3,0)
        colors[p_i] = paste0(colors[p_i], as.hexmode(80+((p_i-1) %% 4)*48))
    }
    pdf(file="Fig5legend.pdf", width=10, height=6)
    vioplot(illus_plots,col = colors,ylim=c(-4,2), xlim=c(1,24))
    
    points(x=1,y=-2, pch=4, col="black", cex=3)
    
    dev.off()
}

true_values = c(0.05333, 0.071862, 0.026892, 1.1107)
plot_figure6 <- function() {
    pdf(file=paste0("Fig_6a.pdf"), width=10, height=6)
    colors = rep(fun_col(length(kfp_subsystems)), each=4)
    for (c_i in seq(length(colors))) {
        colors[c_i] = paste0(colors[c_i], as.hexmode(80+((c_i-1) %% 4)*48))
    }
    vioplot(na.omit(log_kfp_estimate_data),col = colors, ylim=c(-4,1))
    points(rep(log(true_values,10),each=4), pch=4, col="black", cex=3)
    dev.off()
}






calc_input <- function(times, res_sub, met_x, isotopologue_x) {
    text = paste0("calc_function = function(t) { return (",res_sub[["all"]][[met_x]][["anFun"]][[isotopologue_x]][1],")}")
    eval(parse(text=text))
    return(calc_function(times))
}

subs_of_interest = array( c("Gluzzh_1",    "Thrzzc_1",    "Ilezzh_1",    "HxxCyszzc_1",    "Serzzc_1",    "AGNzzh_1",    "Phezzh_1",    "Metzzc_1",    "Aspzzh_1",
                            "Gluzzh_1-M1", "Thrzzc_1-M1", "Ilezzh_1-M1", "HxxCyszzc_1-M1", "Serzzc_1-M1", "AGNzzh_1-M1", "Phezzh_1-M1", "Metzzc_1-M1", "Aspzzh_1-M1"), dim=c(9,2))


.GlobalEnv$use_only_double_logistic <- FALSE
.GlobalEnv$use_only_logistic <- FALSE
.GlobalEnv$use_only_kfp <- FALSE

start_time <- Sys.time()
mc_iter = 20

all_enr_in_logistic = list(list(), list(), list(), list())
all_enr_in_double_logistic = list(list(), list(), list(), list())
all_enr_in_kfp = list(list(), list(), list(), list())

for (mc_run in seq(mc_iter+1)) {
    for (name in names(timepoint_selections)) {
        subset_times = timepoint_selections[[name]]
        print(paste0("Run #",as.character(mc_run), " for  " ,name, " timepoints"))
        enr_input <- array(seq(2*length(subset_times)), dim=c(length(subset_times),2))
        enr_input[,1] = all_enr_in[[name]][[as.character(i)]][["all"]][["Gluzzh_1"]][["exp"]]
        enr_input[,2] = all_enr_in[[name]][[as.character(i)]][["all"]][["Aspzzh_1"]][["exp"]]
        colnames(enr_input) = c("Gluzzh_1", "Aspzzh_1")
        #print(enr_input)
        .GlobalEnv$use_only_double_logistic <- TRUE
        all_enr_in_double_logistic[[name]][[as.character(mc_run)]] <- fit_label_input(enr_input, t=times[subset_times], file="res_fit_enr", mc.cores=NULL)
        .GlobalEnv$use_only_double_logistic <- FALSE
        .GlobalEnv$use_only_kfp <- TRUE
        all_enr_in_kfp[[name]][[as.character(mc_run)]] <- fit_label_input(enr_input, t=times[subset_times], file="res_fit_enr", mc.cores=NULL)
        .GlobalEnv$use_only_kfp <- FALSE
        
    }
}
end_time <- Sys.time()
end_time - start_time



calc_input2 <- function(times, formula) {
    eval(parse(text=paste0("calc_function = function(t) { return (",formula,")}")))
    return(calc_function(times))
}



subs_of_interest = array( c("Gluzzh_1",    "Aspzzh_1",
                            "Gluzzh_1-M1", "Aspzzh_1-M1"), dim=c(2,2))

analyze_alt_an_funs <- function() {
    org_data = res_org$res_dyn$isotopologues
    
    an_dl_fun_res = list()
    an_kfp_fun_res = list()
    for (j in seq(dim(subs_of_interest)[1])) {
        an_fun_diff_sub = subs_of_interest[j,2]
        for (tp in c("6", "11", "21", "81")) {
            #for (tp in c("6")) {
            an_fun_diff_sub_tp = paste0(an_fun_diff_sub,"_",tp)
            for (i in seq(21)) {
                res_current = calc_input2(times, all_enr_in_double_logistic[[tp]][[as.character(i)]][["all"]][[subs_of_interest[j,1]]][["anFun"]][[subs_of_interest[j,2]]] )
                an_dl_fun_res[[an_fun_diff_sub_tp]][[i]] = res_current
                res_current = calc_input2(times, all_enr_in_kfp[[tp]][[as.character(i)]][["all"]][[subs_of_interest[j,1]]][["anFun"]][[subs_of_interest[j,2]]] )
                an_kfp_fun_res[[an_fun_diff_sub_tp]][[i]] = res_current
            }
            pdf(file=paste0("Fig_new_3b_dl_",an_fun_diff_sub_tp,".pdf"), width=10, height=6)
            plot(x=times, y=an_dl_fun_res[[an_fun_diff_sub_tp]][[1]],pch=20, ylim=c(0,1), type="l")
            for (i in seq(2,21)) { lines(x=times, y=an_dl_fun_res[[an_fun_diff_sub_tp]][[i]], pch=20)}
            lines(x=times, y=org_data[,subs_of_interest[j,2]], pch=20, col="green", lwd=4)
            dev.off()
            pdf(file=paste0("Fig_new_3b_kfp_",an_fun_diff_sub_tp,".pdf"), width=10, height=6)
            plot(x=times, y=an_kfp_fun_res[[an_fun_diff_sub_tp]][[1]],pch=20, ylim=c(0,1), type="l")
            for (i in seq(2,21)) { lines(x=times, y=an_kfp_fun_res[[an_fun_diff_sub_tp]][[i]], pch=20)}
            lines(x=times, y=org_data[,subs_of_interest[j,2]], pch=20, col="green", lwd=4)
            dev.off()
            
        }
    }
}


# figure 2C
colors = rep(fun_col(3), each=4)
for (c_i in seq(length(colors))) {
    colors[c_i] = paste0(colors[c_i], as.hexmode(80+((c_i-1) %% 4)*48))
}
pdf(file=paste0("Fig_new_2C.pdf"), width=10, height=6)

vioplot(an_fun_diff_sumsq[c(1,2,3,4,seq(29,36))], ylim=c(0,0.2), col=colors)

dev.off()


pdf(file="Fig_new_2C_Glu_6.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Gluzzh_1-M1_6"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Gluzzh_1-M1_6"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Gluzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Glu_11.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Gluzzh_1-M1_11"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Gluzzh_1-M1_11"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Gluzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Glu_21.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Gluzzh_1-M1_21"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Gluzzh_1-M1_21"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Gluzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Glu_81.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Gluzzh_1-M1_81"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Gluzzh_1-M1_81"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Gluzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Met_6.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Metzzc_1-M1_6"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Metzzc_1-M1_6"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Metzzc_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Met_11.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Metzzc_1-M1_11"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Metzzc_1-M1_11"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Metzzc_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Met_21.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Metzzc_1-M1_21"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Metzzc_1-M1_21"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Metzzc_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Met_81.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Metzzc_1-M1_81"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Metzzc_1-M1_81"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Metzzc_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()


pdf(file="Fig_new_2C_Asp_6.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Aspzzh_1-M1_6"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Aspzzh_1-M1_6"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Aspzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Asp_11.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Aspzzh_1-M1_11"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Aspzzh_1-M1_11"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Aspzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Asp_21.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Aspzzh_1-M1_21"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Aspzzh_1-M1_21"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Aspzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()

pdf(file="Fig_new_2C_Asp_81.pdf", width=10, height=6)
plot(x=times,y=an_fun_res[["Aspzzh_1-M1_81"]][[1]],type="l", ylim=c(0,1))
for (i in seq(2,51)) { lines(x=times,y=an_fun_res[["Aspzzh_1-M1_81"]][[i]], ylim=c(0,1))}
lines(x=times,y=org_data[,"Aspzzh_1-M1"], col="green",lwd=4, ylim=c(0,1))
dev.off()


##### Figure 5

fluxratio_fluxes <- c(v193=0.72611, v242=0.027446)

log_fluxratio_estimate_data = list();
for (fn_tp in names(fluxratio_estimate_data)) { 
    log_fluxratio_estimate_data[[fn_tp]] = log(fluxratio_estimate_data[[fn_tp]],10);  
    log_fluxratio_estimate_data[[fn_tp]][log_fluxratio_estimate_data[[fn_tp]]<(-4)] = -4
}

pdf(file="Fig_5e.pdf", width=10, height=6)
colors = rep(fun_col(2), each=4)
for (c_i in seq(length(colors))) {
    colors[c_i] = paste0(colors[c_i], as.hexmode(80+((c_i-1) %% 4)*48))
}
vioplot(log_fluxratio_estimate_data,col = colors, ylim=c(-4,1))
true_values = rep(log(fluxratio_fluxes,10), each=4)
points(true_values, pch=4, col="black", cex=3)
dev.off()


## Figure 8 signal to noise for r168 Asp AT, Glu -> Met

for (tp in names(timepoint_selections)) {
pdf(file=paste0("Fig_new_8_",tp,".pdf"), width=10, height=6)
plot(x=times,y=org_data[,"Gluzzh_1-M1"],type="l", ylim=c(0,0.5))
for (i in seq(1,51)) { points(x=times[timepoint_selections[[tp]]],y=all_enr_in[[tp]][[as.character(i)]]$all$Gluzzh_1$exp, pch=3, cex=1.5, ylim=c(0,0.5))}
for (i in seq(1,51)) { points(x=times[timepoint_selections[[tp]]],y=all_enr_in[[tp]][[as.character(i)]]$all$Aspzzh_1$exp, pch=3, ylim=c(0,0.5), col="red")}
lines(x=times,y=org_data[,"Aspzzh_1-M1"], col="red", ylim=c(0,0.5))
dev.off()
}


