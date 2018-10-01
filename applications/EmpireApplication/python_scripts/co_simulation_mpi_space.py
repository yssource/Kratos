try:
    import KratosMultiphysics
    from KratosMultiphysics.mpi import mpi
except ImportError:
    raise Exception("Running in MPI currently requires Kratos-MPI!")

class CoSimulationMPISpace(object):
    """This requires some MPI-commands exposed to Python
    This is currently only available with Kratos,
    i.e. MPI can only be used with Kratos compiled with MPI
    """

    def __init__(self):
        # Precompute rank and size such that they don't have to be recomputed all the time
        self.comm_rank = mpi.rank
        self.comm_size = mpi.size

        if self.comm_size < 2:
            raise Exception("Running in MPI requires at least 2 processes!")

    def IsDistributed(self):
        return True

    def Barrier(self):
        mpi.world.barrier()

    def Rank(self):
        return self.comm_rank

    def Size(self):
        return self.comm_size