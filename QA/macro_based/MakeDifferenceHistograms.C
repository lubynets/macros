
// ------------------------------------------------------------
// Process one TTree
// ------------------------------------------------------------

void ProcessDifferenceTree(
    TTree* tree,
    TDirectory* outputDir,
    int nBins = 100,
    double tolerance = 0.0)
{
    const Long64_t nEntries = tree->GetEntries();

    std::cout << "    Tree: " << tree->GetName()
              << " (" << nEntries << " entries)"
              << std::endl;

    TObjArray* branches = tree->GetListOfBranches();

    for (int ib = 0; ib < branches->GetEntriesFast(); ++ib) {

        auto* branch =
            static_cast<TBranch*>(branches->At(ib));

        TObjArray* leaves = branch->GetListOfLeaves();

        for (int il = 0; il < leaves->GetEntriesFast(); ++il) {

            auto* leaf =
                static_cast<TLeaf*>(leaves->At(il));

            const std::string leafName = leaf->GetName();

            // The first macro produces scalar Double_t leaves.
            if (std::string(leaf->GetTypeName()) != "Double_t" ||
                leaf->GetLenStatic() != 1 ||
                leaf->GetLeafCount() != nullptr) {

                throw std::runtime_error(
                    "Expected scalar Double_t leaf: " +
                    std::string(tree->GetName()) + "/" +
                    leafName
                );
            }

            // --------------------------------------------
            // First pass:
            // Determine whether differences are nonzero,
            // and find their minimum and maximum.
            // --------------------------------------------

            bool hasNonZero = false;

            double minValue =
                std::numeric_limits<double>::infinity();

            double maxValue =
                -std::numeric_limits<double>::infinity();

            Long64_t nNonZero = 0;

            for (Long64_t entry = 0; entry < nEntries; ++entry) {

                if (branch->GetEntry(entry) < 0) {
                    throw std::runtime_error(
                        "Error reading branch " +
                        std::string(branch->GetName())
                    );
                }

                const double value = leaf->GetValue();

                if (!std::isfinite(value)) {
                    throw std::runtime_error(
                        "Non-finite value in leaf " + leafName
                    );
                }

                minValue = std::min(minValue, value);
                maxValue = std::max(maxValue, value);

                if (std::abs(value) > tolerance) {
                    hasNonZero = true;
                    ++nNonZero;
                }
            }

            // --------------------------------------------
            // Skip leaves containing only zero differences
            // --------------------------------------------

            if (!hasNonZero) {
                continue;
            }

            // --------------------------------------------
            // Determine histogram range
            // --------------------------------------------

            double xMin = minValue;
            double xMax = maxValue;

            if (xMin == xMax) {

                // Constant nonzero difference.
                double margin = std::max(
                    std::abs(xMin) * 0.1,
                    1e-12
                );

                xMin -= margin;
                xMax += margin;

            } else {

                // Add 5% margin on both sides.
                const double margin = 0.05 * (xMax - xMin);

                xMin -= margin;
                xMax += margin;
            }

            // --------------------------------------------
            // Create histogram
            // --------------------------------------------

            std::string title =
                std::string(tree->GetName()) +
                " : " + leafName +
                ";Difference;Entries";

            auto hist = std::make_unique<TH1D>(
                leafName.c_str(),
                title.c_str(),
                nBins,
                xMin,
                xMax
            );

            // Manage histogram lifetime explicitly.
            hist->SetDirectory(nullptr);

            // --------------------------------------------
            // Second pass:
            // Fill histogram with ALL differences,
            // including zero values.
            // --------------------------------------------

            for (Long64_t entry = 0; entry < nEntries; ++entry) {

                if (branch->GetEntry(entry) < 0) {
                    throw std::runtime_error(
                        "Error reading branch " +
                        std::string(branch->GetName())
                    );
                }

                hist->Fill(leaf->GetValue());
            }

            // --------------------------------------------
            // Save histogram
            // --------------------------------------------

            outputDir->cd();
            hist->Write();

            std::cout
                << "      Saved: " << leafName
                << " | nonzero: " << nNonZero
                << "/" << nEntries
                << " | range: [" << minValue
                << ", " << maxValue << "]"
                << std::endl;
        }
    }
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------

void MakeDifferenceHistograms(
    const char* inputFilename = "differences.root",
    const char* outputFilename = "difference_histograms.root",
    int nBins = 100,
    double tolerance = 0.0)
{
    if (nBins <= 0)
        throw std::runtime_error("nBins must be positive");

    if (!std::isfinite(tolerance) || tolerance < 0.0)
        throw std::runtime_error(
            "tolerance must be finite and non-negative"
        );

    std::unique_ptr<TFile> input(
        TFile::Open(inputFilename, "READ")
    );

    if (!input || input->IsZombie()) {
        throw std::runtime_error(
            std::string("Cannot open input file: ") +
            inputFilename
        );
    }

    std::unique_ptr<TFile> output(
        TFile::Open(outputFilename, "RECREATE")
    );

    if (!output || output->IsZombie()) {
        throw std::runtime_error(
            std::string("Cannot create output file: ") +
            outputFilename
        );
    }

    std::cout << "\nProcessing difference file: "
              << inputFilename << "\n"
              << std::endl;

    // --------------------------------------------------------
    // Iterate over all top-level directories
    // --------------------------------------------------------

    TIter nextDir(input->GetListOfKeys());

    while (auto* dirKey = static_cast<TKey*>(nextDir())) {

        const std::string dirName = dirKey->GetName();

        // Only DF_ directories.
        if (dirName.rfind("DF_", 0) != 0)
            continue;

        TClass* dirClass =
            TClass::GetClass(dirKey->GetClassName());

        if (!dirClass ||
            !dirClass->InheritsFrom(TDirectory::Class()))
            continue;

        TDirectory* inputDir =
            input->GetDirectory(dirName.c_str());

        if (!inputDir)
            throw std::runtime_error(
                "Cannot access directory " + dirName
            );

        std::cout << "Directory: " << dirName << std::endl;

        TDirectory* outputDir =
            output->mkdir(dirName.c_str());

        if (!outputDir)
            throw std::runtime_error(
                "Cannot create directory " + dirName
            );

        // ----------------------------------------------------
        // Iterate over all TTrees in the directory
        // ----------------------------------------------------

        TIter nextTree(inputDir->GetListOfKeys());

        while (auto* treeKey =
                   static_cast<TKey*>(nextTree())) {

            TClass* treeClass =
                TClass::GetClass(treeKey->GetClassName());

            if (!treeClass ||
                !treeClass->InheritsFrom(TTree::Class()))
                continue;

            const std::string treeName = treeKey->GetName();

            auto* tree =
                inputDir->Get<TTree>(treeName.c_str());

            if (!tree)
                throw std::runtime_error(
                    "Cannot access tree " + treeName
                );

            // Create subdirectory named after the tree.
            TDirectory* treeOutputDir =
                outputDir->mkdir(treeName.c_str());

            if (!treeOutputDir)
                throw std::runtime_error(
                    "Cannot create output directory for tree " +
                    treeName
                );

            ProcessDifferenceTree(
                tree,
                treeOutputDir,
                nBins,
                tolerance
            );
        }
    }

    output->Close();

    std::cout
        << "\nFinished. Histograms saved in: "
        << outputFilename
        << std::endl;
}
