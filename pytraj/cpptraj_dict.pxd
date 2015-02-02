    # Traj_PDBfile.h
cdef extern from "Traj_PDBfile.h":
    ctypedef enum PDBWRITEMODE "Traj_PDBfile::PDBWRITEMODE":
        NONEPDBWRITEMODE "Traj_PDBfile::NONE"
        SINGLE "Traj_PDBfile::SINGLE"
        MODEL "Traj_PDBfile::MODEL"
        MULTIPDB "Traj_PDBfile::MULTI"
    # CpptrajFile.h
cdef extern from "CpptrajFile.h":
    ctypedef enum AccessType "CpptrajFile::AccessType":
        READ "CpptrajFile::READ"
        WRITE "CpptrajFile::WRITE"
        APPEND "CpptrajFile::APPEND"
        UPDATE "CpptrajFile::UPDATE"
    # CpptrajFile.h
cdef extern from "CpptrajFile.h":
    ctypedef enum CompressType "CpptrajFile::CompressType":
        NO_COMPRESSION "CpptrajFile::NO_COMPRESSION"
        GZIP "CpptrajFile::GZIP"
        BZIP2 "CpptrajFile::BZIP2"
        ZIP "CpptrajFile::ZIP"
    # CpptrajFile.h
cdef extern from "CpptrajFile.h":
    ctypedef enum FileType "CpptrajFile::FileType":
        UNKNOWN_TYPE "CpptrajFile::UNKNOWN_TYPE"
        STANDARD "CpptrajFile::STANDARD"
        GZIPFILE "CpptrajFile::GZIPFILE"
        BZIP2FILE "CpptrajFile::BZIP2FILE"
        ZIPFILE "CpptrajFile::ZIPFILE"
        MPIFILE "CpptrajFile::MPIFILE"
    # MaskToken.h
cdef extern from "MaskToken.h":
    ctypedef enum MaskTokenType "MaskToken::MaskTokenType":
        OP_NONE "MaskToken::OP_NONE"
        ResNum "MaskToken::ResNum"
        ResName "MaskToken::ResName"
        AtomNum "MaskToken::AtomNum"
        AtomName "MaskToken::AtomName"
        AtomType "MaskToken::AtomType"
        AtomElement "MaskToken::AtomElement"
        SelectAll "MaskToken::SelectAll"
        OP_AND "MaskToken::OP_AND"
        OP_OR "MaskToken::OP_OR"
        OP_NEG "MaskToken::OP_NEG"
        OP_DIST "MaskToken::OP_DIST"
    # Atom.h
cdef extern from "Atom.h":
    ctypedef enum AtomicElementType "Atom::AtomicElementType":
        UNKNOWN_ELEMENT "Atom::UNKNOWN_ELEMENT"
        HYDROGEN "Atom::HYDROGEN"
        BORON "Atom::BORON"
        CARBON "Atom::CARBON"
        NITROGEN "Atom::NITROGEN"
        OXYGEN "Atom::OXYGEN"
        FLUORINE "Atom::FLUORINE"
        PHOSPHORUS "Atom::PHOSPHORUS"
        SULFUR "Atom::SULFUR"
        CHLORINE "Atom::CHLORINE"
        BROMINE "Atom::BROMINE"
        IRON "Atom::IRON"
        CALCIUM "Atom::CALCIUM"
        IODINE "Atom::IODINE"
        MAGNESIUM "Atom::MAGNESIUM"
        COPPER "Atom::COPPER"
        LITHIUM "Atom::LITHIUM"
        POTASSIUM "Atom::POTASSIUM"
        RUBIDIUM "Atom::RUBIDIUM"
        CESIUM "Atom::CESIUM"
        ZINC "Atom::ZINC"
        SODIUM "Atom::SODIUM"
        ALUMINUM "Atom::ALUMINUM"
        ARGON "Atom::ARGON"
        ARSENIC "Atom::ARSENIC"
        SILVER "Atom::SILVER"
        GOLD "Atom::GOLD"
        ASTATINE "Atom::ASTATINE"
        BERYLLIUM "Atom::BERYLLIUM"
        BARIUM "Atom::BARIUM"
        BISMUTH "Atom::BISMUTH"
        CHROMIUM "Atom::CHROMIUM"
        COBALT "Atom::COBALT"
        CADMIUM "Atom::CADMIUM"
        FRANCIUM "Atom::FRANCIUM"
        GALLIUM "Atom::GALLIUM"
        GERMANIUM "Atom::GERMANIUM"
        HELIUM "Atom::HELIUM"
        HAFNIUM "Atom::HAFNIUM"
        MERCURY "Atom::MERCURY"
        INDIUM "Atom::INDIUM"
        IRIDIUM "Atom::IRIDIUM"
        KRYPTON "Atom::KRYPTON"
        MANGANESE "Atom::MANGANESE"
        MOLYBDENUM "Atom::MOLYBDENUM"
        NEON "Atom::NEON"
        NICKEL "Atom::NICKEL"
        NIOBIUM "Atom::NIOBIUM"
        OSMIUM "Atom::OSMIUM"
        PALLADIUM "Atom::PALLADIUM"
        PLATINUM "Atom::PLATINUM"
        LEAD "Atom::LEAD"
        POLONIUM "Atom::POLONIUM"
        RUTHENIUM "Atom::RUTHENIUM"
        RHODIUM "Atom::RHODIUM"
        RHENIUM "Atom::RHENIUM"
        RADON "Atom::RADON"
        RADIUM "Atom::RADIUM"
        SILICON "Atom::SILICON"
        SCANDIUM "Atom::SCANDIUM"
        SELENIUM "Atom::SELENIUM"
        STRONTIUM "Atom::STRONTIUM"
        TIN "Atom::TIN"
        ANTIMONY "Atom::ANTIMONY"
        TITANIUM "Atom::TITANIUM"
        TECHNETIUM "Atom::TECHNETIUM"
        TELLURIUM "Atom::TELLURIUM"
        TANTALUM "Atom::TANTALUM"
        THALLIUM "Atom::THALLIUM"
        VANADIUM "Atom::VANADIUM"
        TUNGSTEN "Atom::TUNGSTEN"
        XENON "Atom::XENON"
        ZIRCONIUM "Atom::ZIRCONIUM"
        YTTRIUM "Atom::YTTRIUM"
        LUTETIUM "Atom::LUTETIUM"
        EXTRAPT "Atom::EXTRAPT"
    # Traj_Mol2File.h
cdef extern from "Traj_Mol2File.h":
    ctypedef enum MOL2WRITEMODE "Traj_Mol2File::MOL2WRITEMODE":
        NONEMOL2WRITEMODE "Traj_Mol2File::NONE"
        # change name to avoid name conflict
        SINGLEMol2File "Traj_Mol2File::SINGLE"
        MOL "Traj_Mol2File::MOL"
        MULTIMOL2 "Traj_Mol2File::MULTI"
    # DataFile.h
cdef extern from "DataFile.h":
    ctypedef enum DataFormatType "DataFile::DataFormatType":
        DATAFILE "DataFile::DATAFILE"
        XMGRACE "DataFile::XMGRACE"
        GNUPLOT "DataFile::GNUPLOT"
        XPLOR "DataFile::XPLOR"
        OPENDX "DataFile::OPENDX"
        REMLOGDATAFILE "DataFile::REMLOG"
        MDOUT "DataFile::MDOUT"
        EVECS "DataFile::EVECS"
        VECTRAJ "DataFile::VECTRAJ"
        UNKNOWN_DATAFORMAT "DataFile::UNKNOWN_DATA"
    # Analysis_Hist.h
cdef extern from "Analysis_Hist.h":
    ctypedef enum NormMode "Analysis_Hist::NormMode":
        NO_NORM "Analysis_Hist::NO_NORM"
        NORM_SUM "Analysis_Hist::NORM_SUM"
        NORM_INT "Analysis_Hist::NORM_INT"
    # Cluster_HierAgglo.h
cdef extern from "Cluster_HierAgglo.h":
    ctypedef enum LINKAGETYPE "Cluster_HierAgglo::LINKAGETYPE":
        SINGLELINK "Cluster_HierAgglo::SINGLELINK"
        AVERAGELINK "Cluster_HierAgglo::AVERAGELINK"
        COMPLETELINK "Cluster_HierAgglo::COMPLETELINK"
    # Mol2File.h
cdef extern from "Mol2File.h":
    ctypedef enum TRIPOSTAG "Mol2File::TRIPOSTAG":
        MOLECULE "Mol2File::MOLECULE"
        ATOMMOL2 "Mol2File::ATOM"
        BOND "Mol2File::BOND"
        SUBSTRUCT "Mol2File::SUBSTRUCT"
    # ReplicaDimArray.h
cdef extern from "ReplicaDimArray.h":
    ctypedef enum RemDimType "ReplicaDimArray::RemDimType":
        UNKNOWNREPDIM "ReplicaDimArray::UNKNOWN"
        TEMPERATURE "ReplicaDimArray::TEMPERATURE"
        PARTIAL "ReplicaDimArray::PARTIAL"
        HAMILTONIAN "ReplicaDimArray::HAMILTONIAN"
        PH "ReplicaDimArray::PH"
    # Box.h
cdef extern from "Box.h":
    ctypedef enum BoxType "Box::BoxType":
        NOBOX "Box::NOBOX"
        ORTHO "Box::ORTHO"
        TRUNCOCT "Box::TRUNCOCT"
        RHOMBIC "Box::RHOMBIC"
        NONORTHO "Box::NONORTHO"
    # TrajectoryFile.h
cdef extern from "TrajectoryFile.h":
    ctypedef enum TrajFormatType "TrajectoryFile::TrajFormatType":
        AMBERNETCDF "TrajectoryFile::AMBERNETCDF"
        AMBERRESTARTNC "TrajectoryFile::AMBERRESTARTNC"
        PDBFILE "TrajectoryFile::PDBFILE"
        MOL2FILE "TrajectoryFile::MOL2FILE"
        CIF "TrajectoryFile::CIF"
        CHARMMDCD "TrajectoryFile::CHARMMDCD"
        GMXTRX "TrajectoryFile::GMXTRX"
        BINPOS "TrajectoryFile::BINPOS"
        AMBERRESTART "TrajectoryFile::AMBERRESTART"
        TINKER "TrajectoryFile::TINKER"
        AMBERTRAJ "TrajectoryFile::AMBERTRAJ"
        SQM "TrajectoryFile::SQM"
        SDF "TrajectoryFile::SDF"
        CONFLIB "TrajectoryFile::CONFLIB"
        UNKNOWN_TRAJ "TrajectoryFile::UNKNOWN_TRAJ"
    # NetcdfFile.h
cdef extern from "NetcdfFile.h":
    ctypedef enum NCTYPE "NetcdfFile::NCTYPE":
        NC_UNKNOWN "NetcdfFile::NC_UNKNOWN"
        NC_AMBERTRAJ "NetcdfFile::NC_AMBERTRAJ"
        NC_AMBERRESTART "NetcdfFile::NC_AMBERRESTART"
cdef extern from "DataSet.h": 
    # DataSet.h
    ctypedef enum DataType "DataSet::DataType":
        UNKNOWN_DATASET "DataSet::UNKNOWN_DATA"
        DOUBLE "DataSet::DOUBLE"
        FLOAT "DataSet::FLOAT"
        INTEGER "DataSet::INTEGER"
        STRING "DataSet::STRING"
        MATRIX_DBL "DataSet::MATRIX_DBL"
        MATRIX_FLT "DataSet::MATRIX_FLT"
        COORDS "DataSet::COORDS"
        VECTOR "DataSet::VECTOR"
        MODES "DataSet::MODES"
        GRID_FLT "DataSet::GRID_FLT"
        REMLOGDATATYPE "DataSet::REMLOG"
        XYMESH "DataSet::XYMESH"
        TRAJ "DataSet::TRAJ"
    # DataSet.h
    ctypedef enum scalarMode "DataSet::scalarMode":
        M_DISTANCE "DataSet::M_DISTANCE"
        M_ANGLE "DataSet::M_ANGLE"
        M_TORSION "DataSet::M_TORSION"
        M_PUCKER "DataSet::M_PUCKER"
        M_RMS "DataSet::M_RMS"
        UNKNOWN_MODE "DataSet::UNKNOWN_MODE"
    # DataSet.h
    ctypedef enum scalarType "DataSet::scalarType":
        ALPHA "DataSet::ALPHA"
        BETA "DataSet::BETA"
        GAMMA "DataSet::GAMMA"
        DELTA "DataSet::DELTA"
        EPSILON "DataSet::EPSILON"
        ZETA "DataSet::ZETA"
        PUCKER "DataSet::PUCKER"
        CHI "DataSet::CHI"
        H1P "DataSet::H1P"
        C2P "DataSet::C2P"
        PHI "DataSet::PHI"
        PSI "DataSet::PSI"
        PCHI "DataSet::PCHI"
        OMEGA "DataSet::OMEGA"
        NOE "DataSet::NOE"
        UNDEFINEDSCALARTYPE "DataSet::UNDEFINED"
    # DataSet_2D.h
cdef extern from "DataSet_2D.h":
    ctypedef enum MatrixType "DataSet_2D::MatrixType":
        NO_OP "DataSet_2D::NO_OP"
        DIST "DataSet_2D::DIST"
        COVAR "DataSet_2D::COVAR"
        MWCOVAR "DataSet_2D::MWCOVAR"
        CORREL "DataSet_2D::CORREL"
        DISTCOVAR "DataSet_2D::DISTCOVAR"
        IDEA "DataSet_2D::IDEA"
        IRED "DataSet_2D::IRED"
        DIHCOVAR "DataSet_2D::DIHCOVAR"
        NMAT "DataSet_2D::NMAT"
    # DataSet_2D.h
cdef extern from "DataSet_2D.h":
    ctypedef enum MatrixKind "DataSet_2D::MatrixKind":
        FULL "DataSet_2D::FULL"
        HALF "DataSet_2D::HALF"
        TRI "DataSet_2D::TRI"
    # ClusterSieve.h
cdef extern from "ClusterSieve.h":
    ctypedef enum SieveType "ClusterSieve::SieveType":
        NONESieveType "ClusterSieve::NONE"
        REGULAR "ClusterSieve::REGULAR"
        RANDOM "ClusterSieve::RANDOM"
    # Dimension.h
cdef extern from "Dimension.h":
    ctypedef enum DimIdxType "Dimension::DimIdxType":
        X "Dimension::X"
        Y "Dimension::Y"
        Z "Dimension::Z"
    # PDBfile.h
cdef extern from "PDBfile.h":
    ctypedef enum PDB_RECTYPE "PDBfile::PDB_RECTYPE":
        ATOMPDB "PDBfile::ATOM"
        HETATM "PDBfile::HETATM"
        CRYST1 "PDBfile::CRYST1"
        TER "PDBfile::TER"
        END "PDBfile::END"
        ANISOU "PDBfile::ANISOU"
        END_OF_FILE "PDBfile::END_OF_FILE"
        UNKNOWNPDBFILE "PDBfile::UNKNOWN"
    # Action.h
cdef extern from "Action.h":
    ctypedef enum RetTypeAct "Action::RetType":
        OKACTION "Action::OK"
        ERRACTION "Action::ERR"
        USEORIGINALFRAME "Action::USEORIGINALFRAME"
        SUPPRESSCOORDOUTPUT "Action::SUPPRESSCOORDOUTPUT"
    # ParmFile.h
cdef extern from "ParmFile.h":
    ctypedef enum ParmFormatType "ParmFile::ParmFormatType":
        AMBERPARM "ParmFile::AMBERPARM"
        PDBFILEPARM "ParmFile::PDBFILE"
        MOL2FILEPARM "ParmFile::MOL2FILE"
        CHARMMPSF "ParmFile::CHARMMPSF"
        CIFFILE "ParmFile::CIFFILE"
        GMXTOP "ParmFile::GMXTOP"
        SDFFILE "ParmFile::SDFFILE"
        # change name to avoid conflict
        TINKERPARM "ParmFile::TINKER"
        UNKNOWN_PARM "ParmFile::UNKNOWN_PARM"
    # TrajinList.h
cdef extern from "TrajinList.h":
    ctypedef enum TrajModeType "TrajinList::TrajModeType":
        UNDEFINEDTRAJINLIST "TrajinList::UNDEFINED"
        NORMAL "TrajinList::NORMAL"
        ENSEMBLE "TrajinList::ENSEMBLE"
    # ClusterList.h
cdef extern from "ClusterList.h":
    ctypedef enum DistModeType "ClusterList::DistModeType":
        USE_FRAMES "ClusterList::USE_FRAMES"
        USE_FILE "ClusterList::USE_FILE"
    # ClusterList.h
cdef extern from "ClusterList.h":
    ctypedef enum DistMetricType "ClusterList::DistMetricType":
        RMS "ClusterList::RMS"
        DME "ClusterList::DME"
        SRMSD "ClusterList::SRMSD"
        DATA "ClusterList::DATA"
    # GridAction.h
cdef extern from "GridAction.h":
    ctypedef enum GridModeType "GridAction::GridModeType":
        ORIGIN "GridAction::ORIGIN"
        BOX "GridAction::BOX"
        MASKCENTER "GridAction::MASKCENTER"
        SPECIFIEDCENTER "GridAction::SPECIFIEDCENTER"
    # Trajin_Multi.h
cdef extern from "Trajin_Multi.h":
    ctypedef enum TargetType "Trajin_Multi::TargetType":
        NONETARGETTYPE "Trajin_Multi::NONE"
        TEMP "Trajin_Multi::TEMP"
        INDICES "Trajin_Multi::INDICES"
        CRDIDX "Trajin_Multi::CRDIDX"
    # Analysis.h
cdef extern from "Analysis.h":
    ctypedef enum RetTypeAna "Analysis::RetType":
        OKANALYSIS "Analysis::OK"
        ERRANALYSIS "Analysis::ERR"
