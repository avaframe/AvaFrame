#!python
# cython: boundscheck=False, wraparound=False, cdivision=True
""" manage Dams in DFA simulation
"""

cpdef (int, double, double, double, double, double, double, double, double, double, double) getWallInteraction(
                                                                          double , double , double,
                                                                          double , double , double,
                                                                          double , double , double,
                                                                          int, double[:], double[:], double[:],
                                                                          double[:], double[:], double[:],
                                                                          double[:], double[:], double[:],
                                                                          int, int, double, int, double, int,
                                                                          double[:,:], double[:,:], double[:,:],
                                                                          double[:,:], double[:,:])

cpdef (int, int, double, double, double, double, double, double, double, double, double) getIntersection(double, double,
                                                                                            double, double,
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            double[:],
                                                                                            int)

cpdef (int, double) linesIntersect(double, double, double, double,
                                   double, double , double , double)
