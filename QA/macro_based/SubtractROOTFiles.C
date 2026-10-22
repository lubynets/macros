// ------------------------------------------------------------
// Error handling
// ------------------------------------------------------------

void Fail(const std::string& message)
{
    throw std::runtime_error(message);
}

// ------------------------------------------------------------
// Find objects of a given ROOT class inside a directory
// ------------------------------------------------------------

std::set<std::string> GetObjectNames(
    TDirectory* dir,
    const TClass* requiredClass,
    const std::string& prefix = "")
{
    std::set<std::string> names;

    TIter next(dir->GetListOfKeys());

    while (auto* key = static_cast<TKey*>(next())) {

        const std::string name = key->GetName();

        if (!prefix.empty() && name.rfind(prefix, 0) != 0)
            continue;

        TClass* cls = TClass::GetClass(key->GetClassName());

        if (cls && cls->InheritsFrom(requiredClass))
            names.insert(name);
    }

    return names;
}

// ------------------------------------------------------------
// Compare sets of names
// ------------------------------------------------------------

void CheckNames(
    const std::set<std::string>& a,
    const std::set<std::string>& b,
    const std::string& context)
{
    if (a == b)
        return;

    std::ostringstream msg;
    msg << "Object name mismatch in " << context << "\n";

    for (const auto& name : a)
        if (!b.count(name))
            msg << "  Missing in second file: " << name << "\n";

    for (const auto& name : b)
        if (!a.count(name))
            msg << "  Missing in first file: " << name << "\n";

    Fail(msg.str());
}

// ------------------------------------------------------------
// Supported numeric types
// ------------------------------------------------------------

bool IsNumericType(const std::string& type)
{
    static const std::set<std::string> types = {
        "Char_t",
        "UChar_t",
        "Short_t",
        "UShort_t",
        "Int_t",
        "UInt_t",
        "Long_t",
        "ULong_t",
        "Long64_t",
        "ULong64_t",
        "Float_t",
        "Double_t",
        "Float16_t",
        "Double32_t"
    };

    return types.count(type) != 0;
}

// ------------------------------------------------------------
// Collect and validate scalar leaves
//
// Each branch must contain exactly one scalar numeric leaf.
// ------------------------------------------------------------

std::map<std::string, TLeaf*> GetLeaves(TTree* tree)
{
    std::map<std::string, TLeaf*> result;

    TObjArray* branches = tree->GetListOfBranches();

    for (int i = 0; i < branches->GetEntriesFast(); ++i) {

        auto* branch =
            static_cast<TBranch*>(branches->At(i));

        TObjArray* leaves = branch->GetListOfLeaves();

        if (leaves->GetEntriesFast() != 1) {
            Fail(
                std::string("Branch ") + branch->GetName() +
                " does not contain exactly one leaf"
            );
        }

        auto* leaf = static_cast<TLeaf*>(leaves->At(0));

        const std::string type = leaf->GetTypeName();

        if (!IsNumericType(type)) {
            Fail(
                std::string("Unsupported leaf type: ") +
                branch->GetName() + " (" + type + ")"
            );
        }

        // Reject variable-length arrays.
        if (leaf->GetLeafCount() != nullptr) {
            Fail(
                std::string("Variable-length array not supported: ") +
                branch->GetName()
            );
        }

        // Reject fixed-size arrays.
        if (leaf->GetLenStatic() != 1) {
            Fail(
                std::string("Array leaf not supported: ") +
                branch->GetName()
            );
        }

        result[branch->GetName()] = leaf;
    }

    return result;
}

// ------------------------------------------------------------
// Validate two trees
// ------------------------------------------------------------

void ValidateTrees(TTree* t1, TTree* t2)
{
    const std::string treeName = t1->GetName();

    if (t1->GetEntries() != t2->GetEntries()) {

        std::ostringstream msg;

        msg << "Different number of entries in tree "
            << treeName
            << ": "
            << t1->GetEntries()
            << " vs "
            << t2->GetEntries();

        Fail(msg.str());
    }

    auto leaves1 = GetLeaves(t1);
    auto leaves2 = GetLeaves(t2);

    std::set<std::string> names1, names2;

    for (const auto& p : leaves1)
        names1.insert(p.first);

    for (const auto& p : leaves2)
        names2.insert(p.first);

    CheckNames(names1, names2, treeName);

    for (const auto& name : names1) {

        TLeaf* l1 = leaves1.at(name);
        TLeaf* l2 = leaves2.at(name);

        if (std::string(l1->GetTypeName()) !=
            std::string(l2->GetTypeName())) {

            Fail(
                "Leaf type mismatch in tree " +
                treeName + ", branch " + name
            );
        }

        if (std::string(l1->GetName()) !=
            std::string(l2->GetName())) {

            Fail(
                "Leaf name mismatch in tree " +
                treeName + ", branch " + name
            );
        }
    }
}

// ------------------------------------------------------------
// Subtract two trees
// ------------------------------------------------------------

void SubtractTrees(
    TTree* t1,
    TTree* t2,
    TDirectory* outputDir)
{
    ValidateTrees(t1, t2);

    auto leaves1 = GetLeaves(t1);
    auto leaves2 = GetLeaves(t2);

    outputDir->cd();

    auto outputTree = std::make_unique<TTree>(
        t1->GetName(),
        t1->GetTitle()
    );

    // Each output branch stores a Double_t.
    //
    // std::map keeps addresses of stored values stable.

    std::map<std::string, Double_t> outputValues;

    for (const auto& item : leaves1) {

        const std::string& branchName = item.first;

        TLeaf* leaf = item.second;

        outputValues[branchName] = 0.0;

        std::string leafList =
            std::string(leaf->GetName()) + "/D";

        outputTree->Branch(
            branchName.c_str(),
            &outputValues.at(branchName),
            leafList.c_str()
        );
    }

    const Long64_t nEntries = t1->GetEntries();

    std::cout
        << "    Processing tree: "
        << t1->GetName()
        << " (" << nEntries << " entries, "
        << leaves1.size() << " branches)"
        << std::endl;

    for (Long64_t entry = 0; entry < nEntries; ++entry) {

        if (t1->GetEntry(entry) < 0 ||
            t2->GetEntry(entry) < 0) {

            Fail(
                "Error reading entry " +
                std::to_string(entry) +
                " in tree " + t1->GetName()
            );
        }

        for (const auto& item : leaves1) {

            const std::string& name = item.first;

            TLeaf* l1 = item.second;
            TLeaf* l2 = leaves2.at(name);

            const Double_t v1 = l1->GetValue();
            const Double_t v2 = l2->GetValue();

            outputValues.at(name) = v1 - v2;
        }

        outputTree->Fill();
    }

    outputTree->Write();
}

// ------------------------------------------------------------
// Main macro
// ------------------------------------------------------------

void SubtractROOTFiles(
    const char* filename1,
    const char* filename2,
    const char* outputFilename = "differences.root")
{
    try {

        // ----------------------------------------------------
        // Open input files
        // ----------------------------------------------------

        std::unique_ptr<TFile> file1(
            TFile::Open(filename1, "READ")
        );

        std::unique_ptr<TFile> file2(
            TFile::Open(filename2, "READ")
        );

        if (!file1 || file1->IsZombie())
            Fail(std::string("Cannot open file: ") + filename1);

        if (!file2 || file2->IsZombie())
            Fail(std::string("Cannot open file: ") + filename2);

        // ----------------------------------------------------
        // Find all DF_ directories
        // ----------------------------------------------------

        auto dirs1 = GetObjectNames(
            file1.get(),
            TDirectory::Class(),
            "DF_"
        );

        auto dirs2 = GetObjectNames(
            file2.get(),
            TDirectory::Class(),
            "DF_"
        );

        CheckNames(dirs1, dirs2, "top-level DF_ directories");

        // ----------------------------------------------------
        // Preflight validation of every directory and tree
        //
        // Do this BEFORE creating the output file.
        // ----------------------------------------------------

        std::cout << "Validating input files..." << std::endl;

        for (const auto& dirName : dirs1) {

            TDirectory* d1 = file1->GetDirectory(dirName.c_str());
            TDirectory* d2 = file2->GetDirectory(dirName.c_str());

            if (!d1 || !d2)
                Fail("Cannot access directory " + dirName);

            auto trees1 = GetObjectNames(
                d1, TTree::Class()
            );

            auto trees2 = GetObjectNames(
                d2, TTree::Class()
            );

            CheckNames(trees1, trees2, dirName);

            for (const auto& treeName : trees1) {

                auto* t1 = d1->Get<TTree>(treeName.c_str());
                auto* t2 = d2->Get<TTree>(treeName.c_str());

                if (!t1 || !t2)
                    Fail("Cannot access tree " + treeName);

                ValidateTrees(t1, t2);
            }
        }

        std::cout << "Validation successful." << std::endl;

        // ----------------------------------------------------
        // Create output file
        // ----------------------------------------------------

        std::unique_ptr<TFile> output(
            TFile::Open(outputFilename, "RECREATE")
        );

        if (!output || output->IsZombie())
            Fail("Cannot create output ROOT file");

        // ----------------------------------------------------
        // Process directories
        // ----------------------------------------------------

        for (const auto& dirName : dirs1) {

            std::cout
                << "Processing directory: "
                << dirName << std::endl;

            TDirectory* d1 = file1->GetDirectory(dirName.c_str());
            TDirectory* d2 = file2->GetDirectory(dirName.c_str());

            TDirectory* outDir = output->mkdir(dirName.c_str());

            if (!outDir)
                Fail("Cannot create output directory " + dirName);

            auto treeNames = GetObjectNames(
                d1, TTree::Class()
            );

            for (const auto& treeName : treeNames) {

                auto* t1 = d1->Get<TTree>(treeName.c_str());
                auto* t2 = d2->Get<TTree>(treeName.c_str());

                SubtractTrees(t1, t2, outDir);
            }
        }

        // ----------------------------------------------------
        // Finish
        // ----------------------------------------------------

        output->Write();
        output->Close();

        std::cout
            << "\nSuccessfully created: "
            << outputFilename
            << std::endl;
    }
    catch (const std::exception& e) {

        std::cerr
            << "\nERROR: "
            << e.what()
            << std::endl;

        throw;
    }
}
