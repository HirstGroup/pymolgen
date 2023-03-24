def relabel_tree1(tree):
    def relabel_node(node, parent_label=float('inf')):
        label = node[0]
        children = node[1:]

        if label < parent_label:
            new_node = [label]
            parent_label = label
        else:
            new_node = []

        sorted_children = sorted(children, key=lambda x: x[0])

        for child in sorted_children:
            new_node.append(relabel_node(child, parent_label))

        return new_node

    return relabel_node(tree)

def relabel_tree2(tree):
    def relabel_node(node, parent_label=float('inf')):
        label = node[0]
        children = node[1:]

        if label < parent_label:
            new_node = [label]
            parent_label = label
        else:
            new_node = []

        sorted_children = sorted(children, key=lambda x: x[0])

        for child in sorted_children:
            new_node.append(relabel_node(child, label))

        return new_node

    return relabel_node(tree)

def relabel_tree3(tree):
    def relabel_helper(node, min_label):
        if isinstance(node, int):
            return min_label, node
        new_children = []
        for child in node:
            new_label, new_child = relabel_helper(child, min_label)
            new_children.append(new_child)
            min_label = new_label + 1
        sorted_children = sorted(new_children, key=lambda x: x[0])
        return sorted_children[0][0], [child[1] for child in sorted_children]
        
    return relabel_helper(tree, 1)[1]

def relabel_helper(node, min_label):
    if isinstance(node, int):
        return min_label, node
    
    new_children = []
    for child in node[1:]:
        new_label, new_child = relabel_helper(child, min_label)
        new_children.append((new_label, new_child))

    sorted_children = sorted(new_children, key=lambda x: (x[0], x[1]))

    new_label = sorted_children[0][0]
    new_children = [x[1] for x in sorted_children]
    return new_label, [node[0]] + new_children


def relabel_tree4(tree):
    return relabel_helper(tree, 1)[1]

def relabel_tree(tree):
    # helper function to recursively relabel the tree
    def relabel_helper(node):
        # if node is a leaf, relabel and return
        if isinstance(node, int):
            return (2, node)

        # otherwise, sort children by label and recursively relabel them
        sorted_children = sorted((relabel_helper(child) for child in node), key=lambda x: (x[0], x[1]))
        new_label = sorted_children[0][0] - 1
        new_children = [child[1] for child in sorted_children]

        return (new_label, new_children)

    # call the helper function on the tree
    new_tree = relabel_helper(tree)[1]
    new_tree[0] = 1  # set the root node to 1
    return new_tree

def depth_first_traversal(tree):
    """
    Perform depth-first traversal of the tree and print the values of the nodes
    """
    # Base case: the tree is empty
    if not tree:
        return
    
    # Print the value of the root node
    print(tree[0], end=' ')
    
    # Recursively traverse the left subtree
    depth_first_traversal(tree[1])
    
    # Recursively traverse the right subtree
    depth_first_traversal(tree[2])


tree = [10, [5], [11, [7, [3], [6]], [4, [2], [1]]]]
newtree = depth_first_traversal(tree)
print(tree)
print(newtree)
