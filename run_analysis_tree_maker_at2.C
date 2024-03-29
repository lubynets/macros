
void run_analysis_tree_maker_at2(TString path     = "",
                                TString dataSet   = "",
                                Int_t nEvents     = 1000,                             
                                TString setupName = "sis100_electron") {
  const std::string system = "Au+Au";  // TODO can we read it automatically?
  const float beam_mom     = 12.;

  // --- Logger settings ----------------------------------------------------
  const TString logLevel     = "INFO";
  const TString logVerbosity = "LOW";
  // ------------------------------------------------------------------------

  // -----   Environment   --------------------------------------------------
  const TString myName = "run_analysis_tree_maker";
  const TString srcDir = gSystem->Getenv("VMCWORKDIR");  // top source directory
  // ------------------------------------------------------------------------

  // -----   In- and output file names   ------------------------------------
  TString traFile = path + dataSet + "/" + dataSet + ".tra.root";
//   TString rawFile = path + dataSet + "/" + dataSet + ".raw.root";
  TString rawFile = path + dataSet + "/" + dataSet + ".event.raw.root";
  TString recFile = path + dataSet + "/" + dataSet + ".rec.root";
  TString parFile = path + dataSet + "/" + dataSet + ".par.root";
  const std::string outFile =
    dataSet.Data() + std::string(".analysistree.root");
    
//   TString unigenFile = "/lustre/cbm/users/ogolosov/mc/generators/dcmqgsm_smm/auau/pbeam12agev/mbias/root/dcmqgsm_" + dataSet + ".root";
  // ------------------------------------------------------------------------

//   // -----   In- and output file names   ------------------------------------
//   TString traFile = "/lustre/cbm/users/lubynets/WORKDIR/cbmrootchain/test.tra.root";
// //   TString rawFile = path + dataSet + "/" + dataSet + ".raw.root";
//   TString rawFile = "/lustre/cbm/users/lubynets/WORKDIR/cbmrootchain/test.raw.root";
//   TString recFile = "/lustre/cbm/users/lubynets/WORKDIR/cbmrootchain/test.reco.root";
//   TString parFile = "/lustre/cbm/users/lubynets/WORKDIR/cbmrootchain/test.par.root";
//   const std::string outFile = "test.analysistree.root";
//   // ------------------------------------------------------------------------

//   // -----   In- and output file names   ------------------------------------
//   TString common_path = "/lustre/cbm/users/ogolosov/mc/cbmsim/master_fr_18.2.1_fs_jun19p1/dcmqgsm_smm/auau/12agev/mbias/sis100_electron/";
//   
//   TString traFile = common_path + "transport/1/1.tra.root";
//   TString rawFile = common_path + "ebe/1/1.raw.root";
//   TString recFile = common_path + "ebe/1/1.reco.root";
//   TString parFile = common_path + "transport/1/1.par.root";
//   const std::string outFile =
//     dataSet.Data() + std::string(".analysistree.root");
//   // ------------------------------------------------------------------------

  // -----   Remove old CTest runtime dependency file  ----------------------
  const TString dataDir  = gSystem->DirName(dataSet);
  const TString dataName = gSystem->BaseName(dataSet);
  const TString testName = ("run_treemaker");
  // ------------------------------------------------------------------------

  // -----   Load the geometry setup   -------------------------------------
  std::cout << std::endl;
  const TString setupFile =
    srcDir + "/geometry/setup/setup_" + setupName + ".C";
  const TString setupFunct = "setup_" + setupName + "()";

  std::cout << "-I- " << myName << ": Loading macro " << setupFile << std::endl;

  gROOT->LoadMacro(setupFile);
  gROOT->ProcessLine(setupFunct);
  CbmSetup* setup = CbmSetup::Instance();

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------

  TString geoTag;
  auto* parFileList = new TList();

  std::cout << "-I- " << myName << ": Using raw file " << rawFile << std::endl;
  std::cout << "-I- " << myName << ": Using parameter file " << parFile
            << std::endl;
  std::cout << "-I- " << myName << ": Using reco file " << recFile << std::endl;

  // -----   Reconstruction run   -------------------------------------------
  auto* run         = new FairRunAna();
  auto* inputSource = new FairFileSource(recFile);
  inputSource->AddFriend(traFile);
  inputSource->AddFriend(rawFile);
  run->SetSource(inputSource);
  run->SetOutputFile(outFile.c_str());
  run->SetGenerateRunInfo(kTRUE);
  // ------------------------------------------------------------------------

  // ----- Mc Data Manager   ------------------------------------------------
  auto* mcManager = new CbmMCDataManager("MCManager", 1);
  mcManager->AddFile(traFile);
  run->AddTask(mcManager);
  // ------------------------------------------------------------------------

  // ---   STS track matching   ----------------------------------------------
  auto* matchTask = new CbmMatchRecoToMC();
  matchTask -> SetIncludeMvdHitsInStsTrack(true);
  run->AddTask(matchTask);
  // ------------------------------------------------------------------------

  auto* KF = new CbmKF();
  run->AddTask(KF);
  // needed for tracks extrapolation
  auto* l1 = new CbmL1("CbmL1", 1, 3);
  if (setup->IsActive(ECbmModuleId::kMvd)) {
    setup->GetGeoTag(ECbmModuleId::kMvd, geoTag);
    const TString mvdMatBudgetFileName =
      srcDir + "/parameters/mvd/mvd_matbudget_" + geoTag + ".root";
    l1->SetMvdMaterialBudgetFileName(mvdMatBudgetFileName.Data());
  }
  if (setup->IsActive(ECbmModuleId::kSts)) {
    setup->GetGeoTag(ECbmModuleId::kSts, geoTag);
    const TString stsMatBudgetFileName =
      srcDir + "/parameters/sts/sts_matbudget_" + geoTag + ".root";
    l1->SetStsMaterialBudgetFileName(stsMatBudgetFileName.Data());
  }
  run->AddTask(l1);

  // --- TRD pid tasks
  if (setup->IsActive(ECbmModuleId::kTrd)) {

    CbmTrdSetTracksPidLike* trdLI =
      new CbmTrdSetTracksPidLike("TRDLikelihood", "TRDLikelihood");
    trdLI->SetUseMCInfo(kTRUE);
    trdLI->SetUseMomDependence(kTRUE);
    run->AddTask(trdLI);
    std::cout << "-I- : Added task " << trdLI->GetName() << std::endl;
    //     ------------------------------------------------------------------------
  }

  // AnalysisTree converter
  auto* man = new CbmConverterManager();
  man->SetSystem(system);
  man->SetBeamMomentum(beam_mom);

  man->SetOutputName(outFile, "rTree");
  //  man->SetOutTreeName("aTree");

  man->AddTask(new CbmSimEventHeaderConverter("SimEventHeader"));
  man->AddTask(new CbmRecEventHeaderConverter("RecEventHeader"));
//   man->AddTask(new CbmSimTracksConverter("SimParticles", "", std::string(unigenFile)));
  man->AddTask(new CbmSimTracksConverter("SimParticles"));

  CbmStsTracksConverter* taskCbmStsTracksConverter =
    new CbmStsTracksConverter("VtxTracks", "SimParticles");
  taskCbmStsTracksConverter->SetIsWriteKFInfo();
  taskCbmStsTracksConverter->SetIsReproduceCbmKFPF();
  man->AddTask(taskCbmStsTracksConverter);

  man->AddTask(new CbmRichRingsConverter("RichRings", "VtxTracks"));
  man->AddTask(new CbmTofHitsConverter("TofHits", "VtxTracks"));
  man->AddTask(new CbmTrdTracksConverter("TrdTracks", "VtxTracks"));
  man->AddTask(new CbmPsdModulesConverter("PsdModules"));

  run->AddTask(man);

  // -----  Parameter database   --------------------------------------------
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  auto* parIo1        = new FairParRootFileIo();
  auto* parIo2        = new FairParAsciiFileIo();
  parIo1->open(parFile.Data());
  parIo2->open(parFileList, "in");
  rtdb->setFirstInput(parIo1);
  rtdb->setSecondInput(parIo2);
  rtdb->setOutput(parIo1);
  rtdb->saveOutput();
  // ------------------------------------------------------------------------

  // -----   Intialise and run   --------------------------------------------
  run->Init();

  std::cout << "Starting run" << std::endl;
  run->Run(0, nEvents);
  // ------------------------------------------------------------------------

  timer.Stop();
  const Double_t rtime = timer.RealTime();
  const Double_t ctime = timer.CpuTime();
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << "Output file is " << outFile << std::endl;
  std::cout << "Parameter file is " << parFile << std::endl;

  printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);

  // -----   CTest resource monitoring   ------------------------------------
  FairSystemInfo sysInfo;
  const Float_t maxMemory = sysInfo.GetMaxMemory();
  std::cout << R"(<DartMeasurement name="MaxMemory" type="numeric/double">)";
  std::cout << maxMemory;
  std::cout << "</DartMeasurement>" << std::endl;
  std::cout << R"(<DartMeasurement name="WallTime" type="numeric/double">)";
  std::cout << rtime;
  std::cout << "</DartMeasurement>" << std::endl;
  const Float_t cpuUsage = ctime / rtime;
  std::cout << R"(<DartMeasurement name="CpuLoad" type="numeric/double">)";
  std::cout << cpuUsage;
  std::cout << "</DartMeasurement>" << std::endl;
  // ------------------------------------------------------------------------

  // -----   Finish   -------------------------------------------------------
  std::cout << " Test passed" << std::endl;
  std::cout << " All ok " << std::endl;
  //   Generate_CTest_Dependency_File(depFile);
  // ------------------------------------------------------------------------

  //  RemoveGeoManager();
}
