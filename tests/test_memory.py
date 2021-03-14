"""
Memory usage tests module
"""
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from optimal_guide_finder import memory_limit

class TestMemoryThresholds(unittest.TestCase):
    """
    Tests for memory set RAM usage thresholds
    """
    def test_below_threshold(self):
        """
        Test that set memory limit does not raise an error when limit is not passed
        """
        mem_size = 0.05 # 50 mb
        test_size = 5 * 1024 * 1024  # 5 mb

        try:
            memory_limit.set_limit(mem_size)
            self._assign_variable_of_size(test_size)
        except MemoryError:
            self.fail("Unexpected MemoryError raised")

    def test_below_threshold_multi_thread(self):
        """
        Test that set memory limit does not raise an error when limit is not passed in multiple threads
        """
        mem_size = 0.05 # 50 mb
        test_size = 5 * 1024 * 1024  # 5 mb

        try:
            memory_limit.set_limit(mem_size)
            with ThreadPoolExecutor(max_workers=1) as executor:
                thread = executor.submit(self._assign_variable_of_size, test_size)
                thread.result(timeout=60)
        except MemoryError:
            self.fail("Unexpected MemoryError raised")

    def test_above_threshold(self):
        """
        Test that set memory limit raises an error when limit is passed
        """
        mem_size = 0.05 # 50 mb
        test_size = 500 * 1024 * 1024  # 500 mb

        with self.assertRaises(MemoryError):
            memory_limit.set_limit(mem_size)
            self._assign_variable_of_size(test_size)

    def test_above_threshold_multi_thread(self):
        """
        Test that set memory limit raises an error when limit is passed in multiple threads
        """
        mem_size = 0.05 # 50 mb
        test_size = 500 * 1024 * 1024  # 500 mb

        with self.assertRaises(MemoryError):
            memory_limit.set_limit(mem_size)
            with ThreadPoolExecutor(max_workers=1) as executor:
                thread = executor.submit(self._assign_variable_of_size, test_size)
                thread.result(timeout=60)

    def _assign_variable_of_size(self, size):
        """
        Assigns a string of the given size

        size {int} --variable size in bytes
        """
        test_var = bytearray(size)
        time.sleep(0.05)
