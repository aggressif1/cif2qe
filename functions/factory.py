"""
Factory pattern implementation for CIF2QE components.
Centralizes object creation and manages dependencies.
"""

from .process_logic.cif_processor import CIFProcessor
from .process_logic.supercell_processor import SupercellProcessor
from .process_logic.supercell_verifier import SupercellVerifier
from .process_logic.data_validator import DataValidator
from .input_logic.miller_index_input import MillerIndexInput
from .process_logic.plane_equation_generator import PlaneEquationGenerator
from .process_logic.atom_plane_assignment import AtomPlaneAssignment
from .process_logic.reference_plane_selector import ReferencePlaneSelector
from .process_logic.plane_comparison_analyzer import PlaneComparisonAnalyzer
from .process_logic.reference_plane_manager import ReferencePlaneManager
from .process_logic.supercell_expander import SupercellExpander
from .process_logic.miller_verification_manager import MillerVerificationManager
from .process_logic.sequence_and_cell_calculator import SequenceAndCellCalculator
from .process_logic.qe_output_generator import QEOutputGenerator
from .process_logic.vacuum_processor import VacuumProcessor
from .process_logic.coordinate_transformer import CoordinateTransformer
from .process_logic.structure_processor import StructureProcessor
from .input_logic.vacuum_layer_input import VacuumLayerInput
from .input_logic.hydrogenation_input import HydrogenationInput
from .process_logic.hydrogen_placer import HydrogenPlacer
from .process_logic.hydrogen_optimizer import HydrogenOptimizer
from .process_logic.writers.qe_writer import QEWriter
from .process_logic.writers.supercell_writer import SupercellWriter
from .process_logic.writers.plane_writer import PlaneWriter
from .output_logic.display_manager_phase1 import Phase1DisplayManager
from .output_logic.display_manager_phase2 import Phase2DisplayManager
from .output_logic.display_manager_phase3 import Phase3DisplayManager
from .output_logic.display_manager_phase4 import Phase4DisplayManager
from .output_logic.display_manager_phase5 import Phase5DisplayManager
from .output_logic.console_writer import ConsoleWriter
from .utils.temporary_file_handler import TemporaryFileHandler


class ComponentFactory:
    """
    Factory class for creating CIF2QE components.
    Centralizes object creation and dependency management.
    """

    @staticmethod
    def create_cif_processor():
        """Create CIF processor with all required dependencies."""
        return CIFProcessor()

    @staticmethod
    def create_supercell_processor():
        """Create supercell processor."""
        return SupercellProcessor()

    @staticmethod
    def create_data_validator():
        """Create data validator."""
        return DataValidator()

    @staticmethod
    def create_qe_writer():
        """Create Quantum ESPRESSO writer."""
        return QEWriter()

    @staticmethod
    def create_console_writer():
        """Create console writer."""
        return ConsoleWriter()

    @staticmethod
    def create_temporary_file_handler():
        """Create temporary file handler."""
        return TemporaryFileHandler()

    @staticmethod
    def create_supercell_verifier():
        """Create supercell verifier."""
        return SupercellVerifier()

    @staticmethod
    def create_supercell_writer():
        """Create supercell writer."""
        return SupercellWriter()

    @staticmethod
    def create_miller_index_input():
        """Create Miller index input handler."""
        return MillerIndexInput()

    @staticmethod
    def create_plane_equation_generator():
        """Create plane equation generator."""
        return PlaneEquationGenerator()

    @staticmethod
    def create_atom_plane_assignment():
        """Create atom plane assignment handler."""
        return AtomPlaneAssignment()

    @staticmethod
    def create_reference_plane_selector():
        """Create reference plane selector."""
        return ReferencePlaneSelector()

    @staticmethod
    def create_plane_comparison_analyzer():
        """Create plane comparison analyzer."""
        return PlaneComparisonAnalyzer()

    @staticmethod
    def create_plane_writer():
        """Create plane writer."""
        return PlaneWriter()

    @staticmethod
    def create_reference_plane_manager():
        """Create Phase3 reference plane manager."""
        return Phase3ReferencePlaneManager()

    @staticmethod
    def create_supercell_expander():
        """Create Phase3 supercell expander."""
        return Phase3SupercellExpander()

    @staticmethod
    def create_miller_verification_manager():
        """Create Phase3 Miller verification manager."""
        return Phase3MillerVerificationManager()

    @staticmethod
    def create_sequence_and_cell_calculator():
        """Create Phase3 sequence and cell calculator."""
        return Phase3SequenceAndCellCalculator()

    @staticmethod
    def create_qe_output_generator():
        """Create QE output generator."""
        return QEOutputGenerator()

    @staticmethod
    def create_phase3_display_manager():
        """Create Phase3 display manager."""
        return Phase3DisplayManager()

    @staticmethod
    def create_vacuum_layer_input():
        """Create vacuum layer input handler."""
        return VacuumLayerInput()

    @staticmethod
    def create_vacuum_processor():
        """Create Phase4 vacuum processor."""
        return Phase4VacuumProcessor()

    @staticmethod
    def create_coordinate_transformer():
        """Create Phase4 coordinate transformer."""
        return Phase4CoordinateTransformer()

    @staticmethod
    def create_structure_processor():
        """Create Phase4 structure processor."""
        return Phase4StructureProcessor()

    @staticmethod
    def create_phase4_display_manager():
        """Create Phase4 display manager."""
        return Phase4DisplayManager()

    @staticmethod
    def create_hydrogenation_input():
        """Create hydrogenation input handler."""
        return HydrogenationInput()

    @staticmethod
    def create_hydrogen_placer():
        """Create Phase5 hydrogen placer."""
        return Phase5HydrogenPlacer()

    @staticmethod
    def create_hydrogen_optimizer():
        """Create Phase5 hydrogen optimizer."""
        return Phase5HydrogenOptimizer()

    @staticmethod
    def create_phase5_display_manager():
        """Create Phase5 display manager."""
        return Phase5DisplayManager()

    @staticmethod
    def create_phase1_display_manager():
        """Create Phase1 display manager."""
        return Phase1DisplayManager()

    @staticmethod
    def create_phase2_display_manager():
        """Create Phase2 display manager."""
        return Phase2DisplayManager()

    @staticmethod
    def create_all_phase1_components():
        """
        Create all components needed for Phase1.
        Returns a dictionary of initialized components.
        """
        return {
            'cif_processor': ComponentFactory.create_cif_processor(),
            'supercell_processor': ComponentFactory.create_supercell_processor(),
            'data_validator': ComponentFactory.create_data_validator(),
            'qe_writer': ComponentFactory.create_qe_writer(),
            'display_manager': ComponentFactory.create_phase1_display_manager(),
            'console_writer': ComponentFactory.create_console_writer(),
            'temp_handler': ComponentFactory.create_temporary_file_handler()
        }

    @staticmethod
    def create_all_phase2_components():
        """
        Create all components needed for Phase2.
        Returns a dictionary of initialized components.
        """
        return {
            'supercell_processor': ComponentFactory.create_supercell_processor(),
            'supercell_verifier': ComponentFactory.create_supercell_verifier(),
            'data_validator': ComponentFactory.create_data_validator(),
            'supercell_writer': ComponentFactory.create_supercell_writer(),
            'display_manager': ComponentFactory.create_phase2_display_manager(),
            'temp_handler': ComponentFactory.create_temporary_file_handler()
        }

    @staticmethod
    def create_all_phase3_components():
        """
        Create all components needed for Phase3.
        Returns a dictionary of initialized components.
        """
        return {
            'miller_input': ComponentFactory.create_miller_index_input(),
            'plane_generator': ComponentFactory.create_plane_equation_generator(),
            'atom_assigner': ComponentFactory.create_atom_plane_assignment(),
            'reference_selector': ComponentFactory.create_reference_plane_selector(),
            'comparison_analyzer': ComponentFactory.create_plane_comparison_analyzer(),
            'plane_writer': ComponentFactory.create_plane_writer(),
            'data_validator': ComponentFactory.create_data_validator(),
            'reference_plane_manager': ComponentFactory.create_reference_plane_manager(),
            'supercell_expander': ComponentFactory.create_supercell_expander(),
            'miller_verification_manager': ComponentFactory.create_miller_verification_manager(),
            'sequence_and_cell_calculator': ComponentFactory.create_sequence_and_cell_calculator(),
            'qe_output_generator': ComponentFactory.create_qe_output_generator(),
            'display_manager': ComponentFactory.create_phase3_display_manager(),
            'temp_handler': ComponentFactory.create_temporary_file_handler()
        }

    @staticmethod
    def create_all_phase4_components():
        """
        Create all components needed for Phase4.
        Returns a dictionary of initialized components.
        """
        return {
            'vacuum_input': ComponentFactory.create_vacuum_layer_input(),
            'vacuum_processor': ComponentFactory.create_vacuum_processor(),
            'coordinate_transformer': ComponentFactory.create_coordinate_transformer(),
            'structure_processor': ComponentFactory.create_structure_processor(),
            'qe_writer': ComponentFactory.create_qe_writer(),
            'display_manager': ComponentFactory.create_phase4_display_manager()
        }

    @staticmethod
    def create_all_phase5_components():
        """
        Create all components needed for Phase5.
        Returns a dictionary of initialized components.
        """
        return {
            'hydrogenation_input': ComponentFactory.create_hydrogenation_input(),
            'hydrogen_placer': ComponentFactory.create_hydrogen_placer(),
            'hydrogen_optimizer': ComponentFactory.create_hydrogen_optimizer(),
            'qe_output_generator': ComponentFactory.create_qe_output_generator(),
            'display_manager': ComponentFactory.create_phase5_display_manager(),
            'temp_handler': ComponentFactory.create_temporary_file_handler()
        }
