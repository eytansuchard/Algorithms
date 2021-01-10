
import numpy as NP


# Stoer - Wagner minimal graph cut algorithm class.
# Reference: https://en.wikipedia.org/wiki/Stoer%E2%80%93Wagner_algorithm
class GRAPH_CUT_CLASS():

    def __init__(self, i_matrix):
        # super(GRAPH_CUT_CLASS, self).__init__() # No need, no parent.

        #  Checking that the graph matrix is two dimensional.
        if len(i_matrix.shape) != 2:
            return

        # Checking that the graph matrix is square.
        if i_matrix.shape[0] != i_matrix.shape[1]:
            return

        # v_dict will keep track of index: set of merged nodes.
        # Note that the indices will change after mergers.
        # This object will be a dictionary of sets of vertices.
        self.v_dict = {}
        # *****************************************************

        # Dictyionary of edges. If (i, k) is a key then (k, i) is a key and they have the same value.
        # Examle { (0,1):0.7, (1,0):0.7, (7,3):1, (3,7):1 }
        self.w_dict = {}
        self.vertices_set = set()

        for i in range(i_matrix.shape[0]):
            # Add a vertex to the set.
            self.vertices_set.add(i)
            for k in range(i):
                if i_matrix[k, i] > 0:
                    # Add edges to the edges dictionary.
                    self.w_dict[(i, k)] = self.w_dict[(k, i)] = i_matrix[k, i]

    # Merge two vertices/nodes in the graph.
    def function_merge(self, m_vertex1, m_vertex2, m_set, m_w_dict):

        '''
        :param m_vertex1: integer idndex of the first vertex to be merged.
        :param m_vertex2: integer idndex of the second vertex to be merged
        :param m_set: the set of vertices before merging.
        :param m_w_dict: dictionay of edges and their weights.
        :return: The set on one split of the cut.
        '''

        # ************************************************************************************
        # Loop through all the vertices that are connected to m_vertex1 or m_vertex2.
        # ************************************************************************************

        for i in m_set:
            if m_vertex1 == i or m_vertex2 == i:
                continue

            # *********************************************************
            # Merge edge (m_vertex1, i) with edge (m_vertex2, i).
            # *********************************************************

            m_sum = m_w_dict.get((m_vertex1, i), 0)
            m_w2 = m_w_dict.get((m_vertex2, i), 0)
            m_sum += m_w2

            # *********************************************************
            # End of Merge edge (m_vertex1, i) with edge (m_vertex2, i)
            # *********************************************************

            if m_sum > 0:
                m_w_dict[(m_vertex1, i)] = m_w_dict[(i, m_vertex1)] = m_sum

            if m_w2 > 0:
                del m_w_dict[(m_vertex2, i)]
                del m_w_dict[(i, m_vertex2)]

        # ************************************************************************************
        # End of 'Loop through all the vertices that are connected to m_vertex1 or m_vertex2'.
        # ************************************************************************************

        # Since the vertices/nodes are merged, delete the edge between them.
        del m_w_dict[(m_vertex1, m_vertex2)]
        del m_w_dict[(m_vertex2, m_vertex1)]

        # Remove the merged vertec.
        m_set.remove(m_vertex2)

        # Each vertex/node can be a combination of previous nodes.
        m_cut_set = self.v_dict[m_vertex1]

        # After the merger, nodes are united.
        self.v_dict[m_vertex1] = self.v_dict[m_vertex1].union(self.v_dict[m_vertex2])

        # No longer keep track of node m_vertex2 so delete node m_vertex2.
        # Note that after recursive mergers, self.v_dict[m_vertex2] may contain several unified nodes.
        del self.v_dict[m_vertex2]

        return m_cut_set

    # Stoer - Wagner minimal graph cut algorithm part - 1.
    def function_get_minimum_phase_cut(self, gmpc_element0, gmpc_set, gmpc_w_dict):
        '''
        :param gmpc_element0: integer. An arbitrary node/vertex index to beign with.
        :param gmpc_set: set of vertices as integer indices.
        :param gmpc_w_dict: { (i,k):weight } key:value dictionary pairs.
        :return:
        gmpc_max_vertex1 = last merged first index.
        gmpc_max_vertex2 = last merged second index,
        gmpc_last_sum =  the cut value to node gmpc_max_vertex2.
        gmpc_set = output set of vertices/nodes after the merger.
        gmpc_w_dict = output { (i,k):weight } key:value dictionary pairs after the merger.
        gmpc_cut_set = the vertices that are left on one side after the cut edges are removed.
        '''

        # In this function each time the vertex with the maximal cut to gmpc_aux_set is added to,
        # gmpc_aux_set_c starting with arbitrary gmpc_element0.
        # The last two vertices in this loop are merged together.

        gmpc_max_vertex1 = -1
        gmpc_max_vertex2 = -1

        if not gmpc_element0 in gmpc_set:
            return None, None, None, None, None, None

        gmpc_aux_set = gmpc_set.copy()
        # Initialize a pool set from which we first look for the vertex most connected to gmpc_element0.
        gmpc_aux_set.remove(gmpc_element0)
        # Initialize the set to be which we weigh the sum of edges.
        gmpc_aux_set_c = set([gmpc_element0])
        # That is used for the last cut between a node in gmpc_aux_set and between the set gmpc_aux_set_c.
        gmpc_last_sum = 0

        while True:

            gmcp_max_vertex = -1
            gmpc_max_sum = 0

            for i in gmpc_aux_set:

                # ************************************************************
                # Calculate the cut between the set gmpc_aux_set_c and node i.
                # ************************************************************

                gmpc_sum = 0

                for k in gmpc_aux_set_c:
                    if (i, k) in gmpc_w_dict.keys():
                        gmpc_sum += gmpc_w_dict[(i, k)]

                # ************************************************************
                # ************************************************************

                # ******************************************************************************
                # Does node i achieves a maximal cut between node i and the set gmpc_aux_set_c ?
                # ******************************************************************************
                if gmpc_max_sum < gmpc_sum:
                    gmpc_max_sum = gmpc_sum
                    gmpc_last_sum = gmpc_sum
                    gmcp_max_vertex = i
                # ******************************************************************************
                # ******************************************************************************

            if gmcp_max_vertex > -1:
                gmpc_max_vertex2 = gmpc_max_vertex1
                gmpc_max_vertex1 = gmcp_max_vertex

            else:
                break

            gmpc_aux_set.remove(gmcp_max_vertex)
            gmpc_aux_set_c.add(gmcp_max_vertex)

        # We have to merge gmpc_max_vertex1 and gmpc_max_vertex2

        gmpc_cut_set = set()

        # *********************************************************************************************************
        # When the loop has more than 2 connected vertices, the elgorithm remembers the last two vertices.
        # *********************************************************************************************************

        if gmpc_max_vertex1 > -1 and gmpc_max_vertex2 > -1:
            if gmpc_last_sum > 0:
                gmpc_cut_set = self.function_merge(gmpc_max_vertex1, gmpc_max_vertex2, gmpc_set, gmpc_w_dict)
        else:
            gmpc_max_vertex2 = gmpc_max_vertex1
            gmpc_max_vertex1 = gmpc_element0
            if gmpc_last_sum > 0:
                gmpc_cut_set = self.function_merge(gmpc_max_vertex1, gmpc_max_vertex2, gmpc_set, gmpc_w_dict)

        # *********************************************************************************************************
        # End of 'When the loop has mode than 2 connected vertices, the elgorithm remembers the last two vertices'.
        # *********************************************************************************************************

        return gmpc_max_vertex1, gmpc_max_vertex2, gmpc_last_sum, gmpc_set, gmpc_w_dict, gmpc_cut_set

    # Stoer - Wagner minimal graph cut algorithm part 2.
    def function_get_minimum_cut(self):
        gmc_set = self.vertices_set.copy()
        gmc_w_dict = self.w_dict.copy()
        gmc_min = None
        gmc_optimal_set = set()

        gmc_n = len(self.vertices_set)
        self.v_dict = {}

        for i in range(gmc_n):
            self.v_dict[i] = set([i])

        while len(gmc_set) > 1:
            gmc_element = gmc_set.pop()
            gmc_set.add(gmc_element)

            # ***********************************************************
            # Look for a last cut in the first loop of Stoer - Wagner algorithm.
            # ***********************************************************

            gmc_vertex1, gmc_vertex2, gmc_sum, gmc_set, gmc_w_dict, gmc_last_cut_set = \
                self.function_get_minimum_phase_cut(gmc_element, gmc_set, gmc_w_dict)

            # Keep track of the last minimal Stoer - Wagner cut.
            if gmc_min is None:
                gmc_min = gmc_sum
                gmc_optimal_set = gmc_last_cut_set

            elif gmc_min > gmc_sum:
                gmc_min = gmc_sum
                gmc_optimal_set = gmc_last_cut_set

            # ***********************************************************
            # End of 'Keep track of the last minimal Stoer - Wagner cut'.
            # ***********************************************************

        print("{:.2f}".format(gmc_min))

        # This is one of the split sets.
        print(gmc_optimal_set)
        print(self.vertices_set - gmc_optimal_set)

def main():

    ma_x = NP.zeros((8, 8), dtype=NP.float)
    ma_x[3, 7] = ma_x[7, 3] = 2
    ma_x[2, 6] = ma_x[6, 2] = 2
    ma_x[3, 2] = ma_x[2, 3] = 4
    ma_x[6, 7] = ma_x[7, 6] = 3
    ma_x[3, 6] = ma_x[6, 3] = 2
    ma_x[1, 2] = ma_x[2, 1] = 3
    ma_x[6, 5] = ma_x[5, 6] = 1
    ma_x[1, 5] = ma_x[5, 1] = 2
    ma_x[1, 0] = ma_x[0, 1] = 2
    ma_x[0, 4] = ma_x[4, 0] = 3
    ma_x[5, 4] = ma_x[4, 5] = 3
    ma_x[1, 4] = ma_x[4, 1] = 2
    ma_x[1, 5] = ma_x[5, 1] = 2

    # Initialize the Stoer - Wagner minimal graph cut algorithm.
    ma_min_graph_cut = GRAPH_CUT_CLASS(ma_x)
    # Call the Stoer - Wagner minimal graph cut algorithm.
    ma_min_graph_cut.function_get_minimum_cut()

    return

if __name__ == '__main__':
    main()
