"""
methods for preconditioning the rate matrix to improve the performance 
of the linear algebra solvers
"""


import numpy as np
import networkx as nx


class MSTSpectralDecomposition(object):
    """spectral decomposition of a graph using the minimum spanning tree in the limit of small T"""
    def __init__(self, Ei, Eij, T=.1):
        self.Ei = Ei
        self.Eij = Eij
        self.T = T
        
        self.sorted_indices = sorted(((E, i) for i, E in self.Ei.items()))
        print(self.sorted_indices)
        self.sorted_indices = [i for E, i in self.sorted_indices]
        
        self.run()

    def compute_mst(self):
        self.graph = nx.Graph()
        for i, E in self.Ei.items():
            self.graph.add_node(i, E=E)
        for (i, j), E in self.Eij.items():
            self.graph.add_edge(i, j, E=E)
        self.mst = nx.minimum_spanning_tree(self.graph, weight='E')
        return self.mst
    
    def get_edge_energy(self, u, v):
        try:
            return self.Eij[(u,v)]
        except KeyError:
            return self.Eij[(v,u)]
    
    def find_u(self, mst, s, u, v):
        u[s] = 0.
        v[s] = 0. 
        for parent, child in nx.bfs_edges(mst, s):
            E = self.get_edge_energy(parent, child)
            if E > u[parent]:
                u[child] = E
            else:
                u[child] = u[parent]
            v[child] = u[child] - self.Ei[child]
    
    def recompute_uv(self, mst, s, u, v):
        unew = dict()
        vnew = dict()
        self.find_u(mst, s, unew, vnew)
        for i, unewi in unew.items():
            if unewi < u[i]:
                u[i] = unewi
            v[i] = min(v[i], u[i] - self.Ei[i])
        u[s] = 0.
        v[s] = 0.
            
    def get_connected_sink(self, mst, s1, u):
        s2 = None
        for parent, child in nx.bfs_edges(mst, s1):
            if u[child] == 0.:
                s2 = child
                break
        if s2 is None:
            print(nx.connected_components(mst))
        assert s2 is not None
        return s2
    
    def find_cutting_edge(self, mst, s1, u):
        u1 = u[s1]
        for parent, child in nx.bfs_edges(mst, s1):
            print("u[parent], u[child]", s1, ":", parent, child, u[parent], u[child], u1)
            if u[child] < u1:
                assert u[parent] == u1
                i, j = child, parent
                break
        E = self.get_edge_energy(i, j)
#        path = nx.shortest_path(mst, s1, s2)
#        print "path", path
#        E, i, j = max( [(self.get_edge_energy(i,j), i, j) for i, j in izip(path, path[1:])] )
        print("cutting edge", i, j, E)
        return E, i, j
    
    def compute_kth_eigenvector_components(self, mst, s):
        """compute the non zero components of the kth eigenvector"""
        evec = dict()
        for i in nx.node_connected_component(mst, s):
            evec[s] = 1.
        return evec
    
    def matricize_eigenvectors(self, evecs):
        nodes = sorted(self.Ei.keys())
        node2index = dict(((n, i) for i, n in enumerate(nodes)))
        eigenvectors = np.zeros([len(nodes), len(nodes)])
        eigenvectors[:,0] = 1.
        for k, evec in evecs.items():
            for node, v in evec.items():
                i = node2index[node]
                eigenvectors[i,k] = v
        return eigenvectors
    
    def run(self):
        mst = self.compute_mst()
        print(self.mst.number_of_edges())
        
        u = dict()
        v = dict()
        evecs = dict()
        delta = np.zeros(len(self.Ei))
        
        k = 0
        s1 = self.sorted_indices[0]
        sinks = set()
        sinks.add(s1)
        u[s1] = 0.
        v[s1] = 0.
        
        self.find_u(mst, s1, u, v)
        for i, ui in u.items():
            print("u v", i, u, v[i])
            
        for k in range(1, len(self.Ei)):
            print("")
            # find the next sink
            s = max(iter(v.keys()), key=lambda i:v[i])
#            sold = slist[-1]
#            print slist
            print("current sink", s)
            print("past sinks", sinks)
            print("# edges", mst.number_of_edges(), mst.number_of_edges(s))
            
            if True:
                s2 = self.get_connected_sink(mst, s, u)
            
            for i, ui in u.items(): print("  u[",i,"] = ", ui)
            
            # find the new cutting edge
            Epq, p, q = self.find_cutting_edge(mst, s, u)
            print("p q Epq", p, q, Epq)
            print("s Es   ", s, self.Ei[s], "sold Esold", s2, self.Ei[s2]) 
            print("eval", "exp(-(", Epq, "-", self.Ei[s], ")/", self.T, ") = ", end=' ') 
            print(np.exp(-(Epq - self.Ei[s]) / self.T))
            
            # save the eigenvalue deltas
            delta[k] = Epq - self.Ei[s]
            
            # remove the cutting edge
            mst.remove_edge(p, q)
            
            # Sk
            evecs[k] = self.compute_kth_eigenvector_components(mst, s)
            print("evecs[",k,"] = ", evecs[k])
            
            # recompute u and v
            self.recompute_uv(mst, s, u, v)
            
            sinks.add(s)
        
        # process eigenvalues
        print(delta)
#        evals = np.cumsum(delta)
        self.eigenvalues = np.exp(-delta / self.T)
        self.eigenvalues[0] = 0.
        print("eigenvalues", self.eigenvalues)
        
        # process eigenvectors
        self.eigenvectors = self.matricize_eigenvectors(evecs)
        print(self.eigenvectors)
        
        





def test1():
    from numpy import exp
    T = 1.
    beta = 1./T
    E3 = 0.;
    E1 = 2.;
    E2 = 3.;
    E12 = 4.1;
    E23 = 5.2;
    k = dict()
    k[(1,2)] = exp(-beta*(E12-E1));
    k[(1,3)] = 0;
    k[(1,1)] = -k[(1,3)] - k[(1,2)];
    k[(2,1)] = exp(-beta*(E12-E2));
    k[(2,3)] = exp(-beta*(E23-E2));
    k[(2,2)] = -k[(2,1)] - k[(2,3)];
    k[(3,1)] = 0;
    k[(3,2)] = exp(-beta*(E23-E3));
    k[(3,3)] = -k[(3,1)] -k[(3,2)];
    
    Ei = dict()
    Ei[1] = E1
    Ei[2] = E2
    Ei[3] = E3
    Eij = dict()
    Eij[(1,2)] = E12
    Eij[(2,3)] = E23
    
    T = 0.1
    spect = MSTSpectralDecomposition(Ei, Eij, T=T)
    get_eigs(Ei, Eij, T=T)


def make_random_energies_complete(nnodes):
    Ei = {}
    Eij = {}
    for i in range(nnodes):
        Ei[i] = np.random.uniform(-1,1)
    for i in range(nnodes):
        for j in range(i):
            Eij[(j,i)] = max(Ei[i], Ei[j]) + np.random.uniform(.1, 1)
    return Ei, Eij 

def make_rate_matrix(Ei, Eij, T=.05):
    n = len(Ei)
    m = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if i == j: continue
            try:
                Ets = Eij[(i,j)]
            except KeyError:
                try:
                    Ets = Eij[(j,i)]
                except KeyError:
                    # there is no edge i,j
                    continue
            m[i,j] = np.exp(-(Ets - Ei[i])/T)
    for i in range(n):
        m[i,i] = - m[i,:].sum()
    return m

def get_eigs(Ei, Eij, T=0.05):
    from pele.utils.hessian import sort_eigs
    m = make_rate_matrix(Ei, Eij, T=T)
    lam, v = np.linalg.eig(m)
    lam, v = sort_eigs(lam, v)
    print("exact eigenvalues", sorted(-lam))
    print("exact eigenvectors")
    print(v)
#    print v

def test2():
    np.random.seed(0)
    Ei, Eij = make_random_energies_complete(4)
    print(list(Ei.items()))
    print(list(Eij.items()))
    T = .02
    spect = MSTSpectralDecomposition(Ei, Eij, T=T)
    
    get_eigs(Ei, Eij, T=T)
    

if __name__ == "__main__":
    test2()
