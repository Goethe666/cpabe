class BinaryTreeNode:
    def __init__(self, node_id, is_leaf=False, user_id=None):
        self.node_id = node_id
        self.is_leaf = is_leaf
        self.user_id = user_id
        self.left = None
        self.right = None
        self.parent = None
        self.polynomial_coeff = None

class CompleteBinaryTree:
    def __init__(self, num_users):
        self.num_users = num_users
        self.nodes = {}
        self.leaves = {}
        self.root = None
        self._build_tree()

    def _build_tree(self):
        if self.num_users == 0:
            return
        leaf_count = 1
        while leaf_count < self.num_users:
            leaf_count *= 2
        total_nodes = 2 * leaf_count - 1
        for i in range(1, total_nodes + 1):
            is_leaf = (i > total_nodes // 2)
            node = BinaryTreeNode(i, is_leaf)
            self.nodes[i] = node
            if i == 1:
                self.root = node
        for i in range(1, total_nodes // 2 + 1):
            left_child_id = 2 * i
            right_child_id = 2 * i + 1
            if left_child_id <= total_nodes:
                self.nodes[i].left = self.nodes[left_child_id]
                self.nodes[left_child_id].parent = self.nodes[i]
            if right_child_id <= total_nodes:
                self.nodes[i].right = self.nodes[right_child_id]
                self.nodes[right_child_id].parent = self.nodes[i]
        leaf_nodes = [node for node in self.nodes.values() if node.is_leaf]
        for i, user_id in enumerate(range(1, self.num_users + 1)):
            if i < len(leaf_nodes):
                leaf_nodes[i].user_id = user_id
                self.leaves[user_id] = leaf_nodes[i]

    def get_path(self, user_id):
        if user_id not in self.leaves:
            return []
        path = []
        current = self.leaves[user_id]
        while current is not None:
            path.append(current)
            current = current.parent
        return path

    def cover_algorithm(self, revoked_users):
        if not revoked_users:
            return [self.root]
        all_users = set(self.leaves.keys())
        non_revoked_users = all_users - revoked_users
        if not non_revoked_users:
            return []
        cover_nodes = []
        covered_users = set()
        def find_cover(node):
            if node is None:
                return
            node_users = self._get_users_under_node(node)
            node_non_revoked = node_users & non_revoked_users
            node_revoked = node_users & revoked_users
            if not node_non_revoked:
                return
            if not node_revoked and node_non_revoked:
                if not (node_non_revoked <= covered_users):
                    cover_nodes.append(node)
                    covered_users.update(node_non_revoked)
                return
            if node.left:
                find_cover(node.left)
            if node.right:
                find_cover(node.right)
        find_cover(self.root)
        return cover_nodes

    def _get_users_under_node(self, node):
        if node.is_leaf:
            return {node.user_id} if node.user_id else set()
        users = set()
        if node.left:
            users.update(self._get_users_under_node(node.left))
        if node.right:
            users.update(self._get_users_under_node(node.right))
        return users