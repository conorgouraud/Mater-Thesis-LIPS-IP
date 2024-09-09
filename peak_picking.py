import numpy as np
import networkx as nx


def build_centurion_tree(list_peaks): #function modified from mass2chem
    '''
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    list_mass_tracks has similar format as list_peaks.
    '''
    d = {}
    for p in list_peaks:
        cent = int(100 * p['mz'])
        if cent in d:
            d[cent].append(p)
        else:
            d[cent] = [p]
    return d


def find_all_matches_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm): # function modified from mass2chem
    '''
    Return matched peaks in mz_centurion_tree by m/z diff within limit_ppm.
    '''
    q = int(query_mz * 100) 

    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-2, q-1, q, q+1, q+2): #changed from +/- 1  to +/- 2 as for higher mass values +1 is less than 5ppm
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            diff = peak['mz'] - query_mz
            if - mz_tol < diff < mz_tol: # not looking for the halogen small mass diffs yet

                results.append(peak)
            
    return results


def find_all_matches_centurion_indexed_list_for_m32(query_mz, mz_centurion_tree, limit_ppm):
    '''
    Return matched peaks in mz_centurion_tree by m/z diff within limit_ppm.
    '''
    q = int(query_mz * 100) 

    second_peak_mz = query_mz - 1.003355

    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-2, q-1, q, q+1, q+2): #changed from +/- 1  to +/- 2 as for higher mass values +1 is less than 5ppm
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            if query_mz < 1000:
                diff = peak['mz'] - query_mz
                if mz_tol < diff < mz_tol: 
                    if peak['mz'] - second_peak_mz > 0.992:
                        results.append(peak)
            elif query_mz > 1000 and abs(peak['mz'] - query_mz) < mz_tol:
                    results.append(peak)
    return results


#to generate lists directy from mzml file, for each rt spectrum
def generate_peak_lists_for_all_rt(exp):
    all_peak_lists = []
    for spectrum_index, spectrum in enumerate(exp):
        peak_list = []
        rt = spectrum.getRT()
        peaks = spectrum.get_peaks()
        for peak_id, (mz, intensity) in enumerate(zip(*peaks)):
            if intensity > 0: 
                peak_dict = {
                    'mz': mz,
                    'rtime': rt,
                    'height': intensity,
                    'id': f"{peak_id}",
                }
                peak_list.append(peak_dict)
        all_peak_lists.append(peak_list)
    return all_peak_lists


def within_tolerance(value1, value2, target_difference, ppm_tolerance):
            # Calculate the ppm difference
            difference = abs(value1 - value2)
            tolerance = ppm_tolerance * value2
            return abs(difference - target_difference) <= tolerance


def get_isotopic_edge_pairs(list_peaks,  # function modified from kiphu
                    mztree, 
                    mz_tolerance_ppm, 
                    search_patterns = [(1.003355, '13C/12C', (0.1, 0.8))], # (0.1, 0.8) is to do with checking relative abundances but this will not be used here
                    check_isotope_ratio = False, # not checking abundances FALSE
                    ):
    '''
    To find all isotope pairs. 

    Input
    =====
    list_peaks: 
        [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'rtime': 654, 'height': 14388.0, 'id': 555}, ...]
    mztree: 
        indexed list_peaks
    mz_tolerance_ppm: 
        ppm tolerance in examining m/z patterns.
    search_patterns: 
        a list in the format of [(mz difference, notion, (ratio low limit, ratio high limit)), ..]
        This can be obtained through search.isotopic_patterns. The ratios are optional, because 
        1) naturally occuring constrains are based on chemical formula;
        2) rules are different when isotope tracers are introduced to the experiments.
        But it's important to have a comprehensive list here for isotope tracing experiments.
    rt_tolerance: 
        tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
        Default intended as 2 seconds.

    Return
    ======
        list of lists of peak pairs that match search_patterns patterns, 
        e.g.[ (195, 206, '13C/12C'), ...]. 
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [  ] 
        for _pair in search_patterns:
            (mass_difference, _) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if check_isotope_ratio and len(_pair) > 2:  # checking abundance ratio
                    (abundance_ratio_min, abundance_ratio_max) = _pair[2]
                    if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                        matched.append( (P1['id'], P2['id'], P1['height'], P2['height']) )
                else:
                    matched.append( (P1['mz'], P2['mz'], P1['height'], P2['height']) )
        signatures += matched
    
    return signatures


def filter_chains(subnetwork_chain_list): #input is a list of lists and each inner list has all the chains from a subnetwork

    #print(subnetwork_chain_list)

    ppm_tolerance = 5e-6
    target_difference = 1.003355
    filtered_chain_groups = []
    for chain_group in subnetwork_chain_list: # accessing groups of specific subnetworks chains
        filtered_chains = []
        for chain in chain_group:
        
            ##conditions that remove chains##
            if len(chain) < 3: # should not happen at this stage
                continue
            
            # Part which gets rid of first peaks which are small
            while len(chain) >= 3 and chain[0][1] < 0.1*chain[1][1]: # if first peak is small remove, do this multiple times if necessary
                chain = chain[1:]
            if len(chain) < 3: # some chains might've shortened to length < 3
                continue

            triple_total_intensity = chain[0][1] + chain[1][1] + chain[2][1] # getting rid of patterns where first three peak intensities sum to <1000
            if triple_total_intensity < 1000:
                continue

            # MAYBE ADD MORE FILTERING

            filtered_chains.append(chain[0:3]) # only taking first 3 peaks

            
        # checking for duplicates in one subnetwork since 4th+ peak has been removed
        if len(filtered_chains) == 1:
            filtered_chain_groups.append(filtered_chains)
        
        if len(filtered_chains) > 1: 
            unique_filtered_chains = []
            seen_keys = set()
            for chain in filtered_chains:
                key = tuple(chain[:3, 0]) 
                if key not in seen_keys:
                    seen_keys.add(key)
                    unique_filtered_chains.append(chain)
            filtered_chain_groups.append(unique_filtered_chains)

    flattened_chain_list = [x for y in filtered_chain_groups for x in y] # flattening list of list of arrays into list of array now that duplicates are removed 

    #print(flattened_chain_list)
            
    return flattened_chain_list


def find_paths_from_sources_to_sinks(G):
        
        sources = [node for node in G.nodes() if G.in_degree(node) == 0]
        sinks = [node for node in G.nodes() if G.out_degree(node) == 0]

        #print(min(sources))
        #print(max(sinks))

        '''if len(G) > 100: 
            print("Skipped a network of:", len(G))

            pos = nx.spring_layout(G)
            nx.draw(G, pos, with_labels=True, node_size=700, arrows=True)
            node_labels = nx.get_node_attributes(G, 'value')
            nx.draw_networkx_labels(G, pos, labels=node_labels)
            plt.show()

            return

            G_undirected = G.to_undirected()
            def check_paths_in_graph(graph):
                # Get all nodes in the graph
                nodes = list(graph.nodes())
                # Generate all possible pairs of nodes
                for node1, node2 in combinations(nodes, 2):
                    if not nx.has_path(graph, node1, node2):
                        print(f"False: No path exists between {node1} and {node2}")
                    else:
                        print(f"True: Path exists between {node1} and {node2}")
               
            check_paths_in_graph(G_undirected)'''

       
        all_paths_with_values = []
        for source in sources:
            for sink in sinks:
                for path in nx.all_simple_paths(G, source=source, target=sink):
                    path_with_values = [[node, G.nodes[node]['value']] for node in path]
                    all_paths_with_values.append(np.array(path_with_values))
        
        return all_paths_with_values



def peaks_to_networks_to_chains(peak_list,  mz_tolerance_ppm, # heavily modified from kiphu
            isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8))]
                    ):
    mztree = build_centurion_tree(peak_list)
    iso_edges = get_isotopic_edge_pairs(peak_list, mztree, mz_tolerance_ppm=mz_tolerance_ppm,
                                        search_patterns=isotope_search_patterns,
                                        check_isotope_ratio=False,
                                        )
    
    edges = []
    G = nx.DiGraph() # directed

    for e in iso_edges:
        G.add_node(e[0], value=e[2])  # m/z and value: intensity
        G.add_node(e[1], value=e[3])

        if e[0] < e[1]: # used for weakly connected graph later
            edges.append((e[0], e[1]))
        else:
            edges.append((e[1], e[0]))

    G.add_edges_from(edges)
    subnetworks = [G.subgraph(c).copy() for c in nx.weakly_connected_components(G)]
    
    subnetwork_chains = [find_paths_from_sources_to_sinks(subnet) for subnet in subnetworks] # subnetwork chains is a list of list of arrays, so its a list of groups of chains, each group corresponds to a subnetwork

    #print(subnetwork_chains)

    # where we look with the larger ppm tolerance for m32 # This is not active right now, look below at return to activate it
    subnetwork_chains_with_more_possible_third_peaks = [ ]
    for chain_group in subnetwork_chains: 
        new_subnetwork_chains = []
        for chain in chain_group:
            #print('chain', chain)
            P2_mass = chain[1][0] #[second peak][mass of peak]
            mass_difference = 1.003355
            tmp = find_all_matches_centurion_indexed_list_for_m32(P2_mass + mass_difference, mztree, mz_tolerance_ppm) 

            for P3 in tmp: # IF tmp IS EMPTY THEN THIS PART IS NOT DONE HENCE ALL CHAINS OF LENGTH < 3 ARE REMOVED
                #print('This is P3:', P3)
                new_chain = chain 
                if len(new_chain) < 3:  # If there's no third row in the chain array
                    new_chain = np.vstack([new_chain, [P3['mz'], P3['height']]])  # Append a new row with P3 values
                else:
                    new_chain[2] = [P3['mz'], P3['height']]  # Update the existing third peak with P3 values

                #print('new chain', new_chain)
                new_subnetwork_chains.append(new_chain)

        if new_subnetwork_chains: 
            subnetwork_chains_with_more_possible_third_peaks.append(new_subnetwork_chains) 

    
    #filtered_subnetwork_chains = filter_chains(subnetwork_chains_with_more_possible_third_peaks) ## IF YOU WANT TO INCLUDE THE M32 EXTRA PART
    filtered_subnetwork_chains = filter_chains(subnetwork_chains)

    '''for subnetwork in subnetworks: # just to visualise the subnetworks
        if len(subnetwork) > 4:
            print('CHAINS', find_paths_from_sources_to_sinks(subnetwork))
            pos = nx.spring_layout(subnetwork)
            nx.draw(subnetwork, pos, with_labels=True, node_size=700, arrows=True)
            node_labels = nx.get_node_attributes(subnetwork, 'value')
            nx.draw_networkx_labels(subnetwork, pos, labels=node_labels)
            plt.show()'''

    #print(subnetwork_chains)
    #print(filtered_subnetwork_chains)
    return filtered_subnetwork_chains






def is_in_peak(rt, mz, peaks): # function to check if a given rt, mz is in a Peak class
    rt = float(rt)
    mz = float(mz)
    for peak in peaks:
        if peak.rt_min <= rt <= peak.rt_max and peak.mz_min <= mz <= peak.mz_max:
            return peak.id, peak.mz, peak.rt_min, peak.rt_max, peak.rt  # return the centwave EIC ID of the matching peak, the peak mz value, rt_min, and rt_max
    return None, None, None, None, None


# used to make looking up centwave peaks faster
class Peak: 
    def __init__(self, rt, rt_min, rt_max, mz_min, mz_max, mz, id):
        self.rt = float(rt)
        self.rt_min = float(rt_min)
        self.rt_max = float(rt_max)
        self.mz_min = float(mz_min)
        self.mz_max = float(mz_max)
        self.mz = float(mz)
        self.id = id  # keeping EIC id made by centWave

    def __repr__(self):
        return f"Peak(rt={self.rt}, rt_min={self.rt_min}, rt_max={self.rt_max}, mz_min={self.mz_min}, mz_max={self.mz_max}, mz={self.mz}, id={self.id})"





