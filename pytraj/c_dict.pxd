# distutils: language = c++

'''keywords in cpptraj
'''
cdef extern from "Traj_PDBfile.h":
    ctypedef enum PDBWRITEMODE "Traj_PDBfile::PDBWRITEMODE":
        NONE "Traj_PDBfile::NONE"
        SINGLE "Traj_PDBfile::SINGLE"
        MODEL "Traj_PDBfile::MODEL"
        MULTI "Traj_PDBfile::MULTI"


cdef extern from "MetaData.h":
    ctypedef enum scalarType "MetaData::scalarType":
        ALPHA "MetaData::ALPHA"
        BETA "MetaData::BETA"
        GAMMA "MetaData::GAMMA"
        DELTA "MetaData::DELTA"
        EPSILON "MetaData::EPSILON"
        ZETA "MetaData::ZETA"
        NU1 "MetaData::NU1"
        NU2 "MetaData::NU2"
        H1P "MetaData::H1P"
        C2P "MetaData::C2P"
        CHIN "MetaData::CHIN"
        PHI "MetaData::PHI"
        PSI "MetaData::PSI"
        CHIP "MetaData::CHIP"
        OMEGA "MetaData::OMEGA"
        PUCKER "MetaData::PUCKER"
        NOE "MetaData::NOE"
        DIST "MetaData::DIST"
        COVAR "MetaData::COVAR"
        MWCOVAR "MetaData::MWCOVAR"
        CORREL "MetaData::CORREL"
        DISTCOVAR "MetaData::DISTCOVAR"
        IDEA "MetaData::IDEA"
        IREDMAT "MetaData::IREDMAT"
        DIHCOVAR "MetaData::DIHCOVAR"
        IREDVEC "MetaData::IREDVEC"
        UNDEFINED "MetaData::UNDEFINED"


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
        NONE "Traj_Mol2File::NONE"
        # change name to avoid name conflict
        SINGLEMol2File "Traj_Mol2File::SINGLE"
        MOL "Traj_Mol2File::MOL"
        MULTI "Traj_Mol2File::MULTI"

cdef extern from "Mol2File.h":
    ctypedef enum TRIPOSTAG "Mol2File::TRIPOSTAG":
        MOLECULE "Mol2File::MOLECULE"
        ATOMMOL2 "Mol2File::ATOM"
        BOND "Mol2File::BOND"
        SUBSTRUCT "Mol2File::SUBSTRUCT"

cdef extern from "Box.h":
    ctypedef enum BoxType "Box::BoxType":
        NOBOX "Box::NOBOX"
        ORTHO "Box::ORTHO"
        TRUNCOCT "Box::TRUNCOCT"
        RHOMBIC "Box::RHOMBIC"
        NONORTHO "Box::NONORTHO"

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
        GRID_DBL "DataSet::GRID_DBL"
        REMLOGDATATYPE "DataSet::REMLOG"
        XYMESH "DataSet::XYMESH"
        TRAJ "DataSet::TRAJ"
        REF_FRAME "DataSet::REF_FRAME"
        MAT3X3 "DataSet::MAT3X3"
        TOPOLOGY "DataSet::TOPOLOGY"

cdef extern from "Action.h":
    ctypedef enum RetTypeAct "Action::RetType":
        OKACTION "Action::OK"
        ERRACTION "Action::ERR"
        USEORIGINALFRAME "Action::USEORIGINALFRAME"
        SUPPRESSCOORDOUTPUT "Action::SUPPRESSCOORDOUTPUT"

cdef extern from "Analysis.h":
    ctypedef enum RetTypeAna "Analysis::RetType":
        OKANALYSIS "Analysis::OK"
        ERRANALYSIS "Analysis::ERR"
