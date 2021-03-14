"""
Module for setting memory thresholds cross-platform

Based on stackoverflow thread:
https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
"""
import os
import ctypes
import platform

PROCESS_SET_QUOTA = 0x100
PROCESS_TERMINATE = 0x1
JOB_OBJECT_LIMIT_PROCESS_MEMORY = 0x100
JOB_OBJECT_EXTENDED_LIMIT_INFORMATION = 9

class IO_COUNTERS(ctypes.Structure):
    _fields_ = [('ReadOperationCount', ctypes.c_uint64),
                ('WriteOperationCount', ctypes.c_uint64),
                ('OtherOperationCount', ctypes.c_uint64),
                ('ReadTransferCount', ctypes.c_uint64),
                ('WriteTransferCount', ctypes.c_uint64),
                ('OtherTransferCount', ctypes.c_uint64)]

class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [('PerProcessUserTimeLimit', ctypes.c_int64),
                ('PerJobUserTimeLimit', ctypes.c_int64),
                ('LimitFlags', ctypes.c_uint32),
                ('MinimumWorkingSetSize', ctypes.c_void_p),
                ('MaximumWorkingSetSize', ctypes.c_void_p),
                ('ActiveProcessLimit', ctypes.c_uint32),
                ('Affinity', ctypes.c_void_p),
                ('PriorityClass', ctypes.c_uint32),
                ('SchedulingClass', ctypes.c_uint32)]

class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [('BasicLimitInformation', JOBOBJECT_BASIC_LIMIT_INFORMATION),
                ('IoInfo', IO_COUNTERS),
                ('ProcessMemoryLimit', ctypes.c_void_p),
                ('JobMemoryLimit', ctypes.c_void_p),
                ('PeakProcessMemoryUsed', ctypes.c_void_p),
                ('PeakJobMemoryUsed', ctypes.c_void_p)]

def set_limit_windows(pid, size):
    """
    Set process memory limit on windowss

    Arguments:
        pid {int} -- Process ID
        size {int} -- Memory limit in bytes
    """
    job_info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
    out_size = ctypes.c_uint32()

    job = ctypes.windll.kernel32.CreateJobObjectW(None, None)
    assert job != 0

    success = ctypes.windll.kernel32.QueryInformationJobObject(job,
                                                               JOB_OBJECT_EXTENDED_LIMIT_INFORMATION,
                                                               ctypes.POINTER(JOBOBJECT_EXTENDED_LIMIT_INFORMATION)(job_info),
                                                               ctypes.sizeof(JOBOBJECT_EXTENDED_LIMIT_INFORMATION),
                                                               ctypes.POINTER(ctypes.c_uint32)(out_size))
    assert success

    job_info.BasicLimitInformation.LimitFlags |= JOB_OBJECT_LIMIT_PROCESS_MEMORY
    job_info.ProcessMemoryLimit = size
    success = ctypes.windll.kernel32.SetInformationJobObject(job,
                                                             JOB_OBJECT_EXTENDED_LIMIT_INFORMATION,
                                                             ctypes.POINTER(JOBOBJECT_EXTENDED_LIMIT_INFORMATION)(job_info),
                                                             ctypes.sizeof(JOBOBJECT_EXTENDED_LIMIT_INFORMATION))
    assert success

    process = ctypes.windll.kernel32.OpenProcess(PROCESS_SET_QUOTA | PROCESS_TERMINATE,
                                                 False, pid)
    assert process != 0

    success = ctypes.windll.kernel32.AssignProcessToJobObject(job, process)
    assert success


    success = ctypes.windll.kernel32.CloseHandle(job)
    assert success

    success = ctypes.windll.kernel32.CloseHandle(process)
    assert success

def set_limit_unix(size):
    """
    Set process memory limit on windows

    Arguments:
        size {int} -- Memory limit in bytes
    """
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (size, hard))

def set_limit(size):
    """
    Set process memory limit and throw a MemoryError when exceeded

        size {int} -- Memory limit in GiB
    """
    size_bytes = int(size * 1024 * 1024 * 1024)
    if platform.system() == "Windows":
        set_limit_windows(os.getpid(), size_bytes)
    else:
        set_limit_unix(size_bytes)
