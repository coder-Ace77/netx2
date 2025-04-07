import networkx as nx
import random
import os
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict

def generate_graph():
    while True:
        n = random.randint(10,50)
        m = random.randint(1, min(20, n-1))
        G = nx.barabasi_albert_graph(n, m)
        if nx.is_connected(G):
            return G

def stress_centrality_brute_force(G):
    print("BRUTE RUNNING!!")
    """Calculate stress centrality using brute-force approach"""
    stress_centrality = defaultdict(int)
    nodes = sorted(G.nodes())  
    
    for i, source in enumerate(nodes):
        for target in nodes[i+1:]:  
            try:
                paths = list(nx.all_shortest_paths(G, source=source, target=target))
                for path in paths:
                    for node in path[1:-1]:  # Exclude source and target
                        stress_centrality[node] += 1
            except nx.NetworkXNoPath:
                continue
    print("BRUTE STOPPED!!")
    return stress_centrality

def run_test_case(G, test_num,compile_result):
    """Run a single test case and verify results"""
    # Save graph to file
    nx.write_adjlist(G, "grapht.adjlist")
    
    if compile_result.returncode != 0:
        print(f"Test {test_num}: Compilation failed!")
        print(compile_result.stderr)
        return False
    
    run_result = subprocess.run(["./main_exec"], 
                               capture_output=True, text=True)
    if run_result.returncode != 0:
        print(f"Test {test_num}: C++ program execution failed!")
        print(run_result.stderr)
        return False
    
    reference_result = stress_centrality_brute_force(G)
    res = ""
    nodes = sorted(G.nodes())
    for node in nodes:
        res += f"{reference_result[node]} "
    cpp_result = run_result.stdout
    
    passed = True
    if cpp_result!=res:
        passed = False
    
    if passed:
        print(f"Test {test_num}: PASSED")
    else:
        print(f"Test {test_num}: FAILED")
        print("C++ program output:")
        print(run_result.stdout)
    return passed

def main():
    num_tests = 100
    passed_count = 0

    compile_result = subprocess.run(["g++", "-fopenmp", "main.cpp", "-o", "main_exec"], 
                                  capture_output=True, text=True)
    
    for i in range(num_tests):
        G = generate_graph()
        if run_test_case(G, i+1,compile_result):
            passed_count += 1
    
    print(f"\nTest Summary: {passed_count}/{num_tests} passed")
    
    # Clean up
    if os.path.exists("grapht.adjlist"):
        os.remove("grapht.adjlist")
    if os.path.exists("main_exec"):
        os.remove("main_exec")

if __name__ == "__main__":
    main()