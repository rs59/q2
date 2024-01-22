import random

def generate_test_metis_file_verbose(X, Y, M, Q):

    # Generate a random graph as a list of edges
    all_edges = []
    for _ in range(Y):
        start_vertex = random.randint(1, X)
        end_vertex = random.randint(1, X)
        if _%1000==0:
          print(f"{_} of {Y} generated")
        while start_vertex == end_vertex or (start_vertex, end_vertex) in all_edges or (end_vertex, start_vertex) in all_edges:
            end_vertex = random.randint(1, X)
        all_edges.append((start_vertex, end_vertex))
    print(f"Random graph generated of length {len(all_edges)}")

    # Shuffle the edges to randomize them
    random.shuffle(all_edges)
    print(f"Random graph shuffled")

    # Create a disjoint-set data structure to check connectivity
    parent = list(range(X + 1))

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            parent[root_x] = root_y

    # Generate a random spanning tree using Kruskal's algorithm
    spanning_tree_edges = []
    for edge in all_edges:
        start_vertex, end_vertex = edge
        if find(start_vertex) != find(end_vertex):
            union(start_vertex, end_vertex)
            spanning_tree_edges.append(edge)
    print(f"Random spanning tree generated")

    # Generate additional random edges while preserving connectivity
    extra_edges = Y - len(spanning_tree_edges)
    for _ in range(extra_edges):
        start_vertex = random.randint(1, X)
        end_vertex = random.randint(1, X)
        if _%1000==0:
          print(f"{_} of {extra_edges} generated")
        while start_vertex == end_vertex or (start_vertex, end_vertex) in spanning_tree_edges or (end_vertex, start_vertex) in spanning_tree_edges:
            end_vertex = random.randint(1, X)
        spanning_tree_edges.append((start_vertex, end_vertex))
    print(f"Additional random edges generated")

    # Generate all_edge_pairs with weights
    all_edge_pairs = []
    for edge in spanning_tree_edges:
        start_vertex, end_vertex = edge
        edge_weight = random.randint(1, Q)
        all_edge_pairs.append((start_vertex, end_vertex, edge_weight))

    with open(f"x{X}y{Y}m{M}q{Q}.metis", "w") as file:
        # Write the first line with X, Y, and Z
        file.write(f"{X} {Y} 011 1\n")
        print(f"Writing header line: {X} {Y} 3 1") #3 (011 flags) 1: weighted edges and 1 weight per vertex

        for vertex in range(1, X + 1):
            # Generate a random vertex weight between 1 and M
            vertex_weight = random.randint(1, M)
            print(f"Generating vertex {vertex}: Weight={vertex_weight}")

            # Get BN CN pairs for this vertex without modifying all_edge_pairs
            vertex_edge_pairs = []
            for start_vertex, end_vertex, edge_weight in all_edge_pairs:
                if start_vertex == vertex:
                    vertex_edge_pairs.extend([end_vertex, edge_weight])
                elif end_vertex == vertex:
                    vertex_edge_pairs.extend([start_vertex, edge_weight])

            vertex_edge_pairs_str = ' '.join(map(str, vertex_edge_pairs))
            file_line = f"{vertex_weight} {vertex_edge_pairs_str}\n"
            file.write(file_line)
            print(f"Writing vertex {vertex}: {vertex_weight} {vertex_edge_pairs_str}")

if __name__ == "__main__":
    # Get inputs from the user
    print("METIS-format graph generator")
    print("X = # of vertices")
    print("Y = # of edges")
    print("M = upper bound for vertex weight")
    print("Q = upper bound for edge weight")
    X, Y, M, Q = map(int, input("Enter X Y M Q: ").split())

    print("Generating test.metis file...")
    # Generate the test.metis file
    generate_test_metis_file_verbose(X, Y, M, Q)
    print("test.metis file has been generated.")
