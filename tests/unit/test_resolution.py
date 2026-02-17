"""Unit tests for pathSQE.Resolution module"""

import pytest
import numpy as np

from pathSQE.Resolution import Instrument, Resolution


class TestInstrument:
    """Tests for the Instrument class"""

    def test_instrument_initialization_cncs(self):
        """Test CNCS instrument initialization"""
        inst = Instrument("CNCS", 0.05, 3.0)
        assert inst.instName == ("CNCS",)
        assert inst.FWHM == 0.05
        assert inst.ei == 3.0
        assert inst.type == "INS"
        assert inst.L1 == 28.3
        assert inst.L2 == 1.487
        assert inst.L3 == 3.5

    def test_instrument_initialization_arcs(self):
        """Test ARCS instrument initialization"""
        inst = Instrument("ARCS", 0.02, 25.0)
        assert inst.FWHM == 0.02
        assert inst.ei == 25.0
        assert inst.type == "INS"
        assert inst.L1 == 11.6
        assert inst.L2 == 1.99
        assert inst.L3 == 3.0

    def test_instrument_initialization_seq(self):
        """Test SEQ instrument initialization"""
        inst = Instrument("SEQ", 0.1, 14.0)
        assert inst.FWHM == 0.1
        assert inst.ei == 14.0
        assert inst.type == "INS"
        assert inst.L1 == 18.0
        assert inst.L2 == 2.0
        assert inst.L3 == 5.5

    def test_instrument_initialization_herix(self):
        """Test HERIX instrument initialization"""
        inst = Instrument("HERIX", 1.5, 17978.0)
        assert inst.FWHM == 1.5
        assert inst.ei == 17978.0
        assert inst.type == "IXS"

    def test_set_pulse_width(self):
        """Test setPulseWidth method"""
        inst = Instrument("ARCS", 0.02, 25.0)
        original_dtp = inst.dtp
        new_dtp = 0.0000025
        inst.setPulseWidth(new_dtp)
        assert inst.dtp == new_dtp
        assert inst.dtp != original_dtp

    def test_set_chopper_uncertainty(self):
        """Test setChopperUncertainty method"""
        inst = Instrument("ARCS", 0.02, 25.0)
        original_dtc = inst.dtc
        new_dtc = 0.000005
        inst.setChopperUncertainty(new_dtc)
        assert inst.dtc == new_dtc
        assert inst.dtc != original_dtc

    def test_set_sample_uncertainty(self):
        """Test setSampleUncertainty method"""
        inst = Instrument("ARCS", 0.02, 25.0)
        original_dtd = inst.dtd
        new_dtd = 0.000005
        inst.setSampleUncertainty(new_dtd)
        assert inst.dtd == new_dtd
        assert inst.dtd != original_dtd

    def test_instrument_defaults_preserved(self):
        """Test that instrument defaults are correctly set"""
        inst_arcs = Instrument("ARCS", 0.02, 25.0)
        inst_cncs = Instrument("CNCS", 0.05, 3.0)

        # ARCS defaults
        assert inst_arcs.dtp == 0.0000012
        assert inst_arcs.dtd == 0.0000024
        assert inst_arcs.dtc == 0.0000036

        # CNCS defaults
        assert inst_cncs.dtp == 0.00000103
        assert inst_cncs.dtd == 0.000014
        assert inst_cncs.dtc == 0.00000026

    def test_instrument_different_fwhm_values(self):
        """Test instrument with different FWHM values"""
        fwhm_values = [0.01, 0.05, 0.1, 1.5]
        for fwhm in fwhm_values:
            inst = Instrument("ARCS", fwhm, 25.0)
            assert inst.FWHM == fwhm

    def test_instrument_different_incident_energies(self):
        """Test instrument with different incident energies"""
        ei_values = [3.0, 8.0, 14.0, 25.0, 100.0]
        for ei in ei_values:
            inst = Instrument("ARCS", 0.02, ei)
            assert inst.ei == ei


class TestResolution:
    """Tests for the Resolution class"""

    def test_resolution_initialization_constant(self):
        """Test Resolution initialization with constant type"""
        inst = Instrument("ARCS", 0.02, 25.0)
        sigma_q = 0.01
        res = Resolution("constant", inst, sigmaq=sigma_q)

        assert res._name == "constant"
        assert res.inst == inst
        assert res.sigmaq == sigma_q

    def test_resolution_initialization_polynomial(self):
        """Test Resolution initialization with polynomial type"""
        inst = Instrument("ARCS", 0.02, 25.0)
        poly_coeffs = [1.0, 0.5, 0.1, 0.01, 0.001]
        res = Resolution("polynomial", inst, poly=poly_coeffs)

        assert res._name == "polynomial"
        assert res.inst == inst
        assert res.poly == poly_coeffs

    def test_resolution_with_all_parameters(self):
        """Test Resolution initialization with all parameters"""
        inst = Instrument("CNCS", 0.05, 3.0)
        sigma_q = 0.02
        poly_coeffs = [1.0, 0.2, 0.05]
        res = Resolution("hybrid", inst, sigmaq=sigma_q, poly=poly_coeffs)

        assert res.sigmaq == sigma_q
        assert res.poly == poly_coeffs

    def test_resolution_without_optional_params(self):
        """Test Resolution initialization without optional parameters"""
        inst = Instrument("ARCS", 0.02, 25.0)
        res = Resolution("constant", inst)

        assert res.sigmaq is None
        assert res.poly is None

    def test_resolution_different_instruments(self):
        """Test Resolution with different instruments"""
        instruments = ["ARCS", "CNCS", "SEQ", "HERIX"]
        for inst_name in instruments:
            inst = Instrument(inst_name, 0.1, 10.0)
            res = Resolution("constant", inst, sigmaq=0.01)
            assert res.inst.instName == (inst_name,)

    def test_resolution_instrument_reference(self):
        """Test that Resolution maintains reference to instrument"""
        inst = Instrument("ARCS", 0.02, 25.0)
        res = Resolution("constant", inst, sigmaq=0.01)

        # Modify instrument parameter
        inst.setPulseWidth(0.0000015)

        # Check that resolution's instrument reference is updated
        assert res.inst.dtp == 0.0000015
