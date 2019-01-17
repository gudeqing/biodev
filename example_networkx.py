from matplotlib import pyplot as plt
import networkx as nx


class StateGraph(object):
    def __init__(self, state):
        """
        drawing
        :param state: state dict from RunCommands.state
        """
        self.state = state
        self.graph = nx.DiGraph()

    def add_edges(self):
        for target in self.state:
            sources = self.state[target]['depend'].strip()
            if sources:
                sources = sources.split(',')
                edges = zip(sources, [target]*len(sources))
                self.graph.add_edges_from(edges)
            else:
                self.graph.add_edge('Input', target, color='green')

    def get_color_dict(self):
        colors = list()
        all_nodes = self.graph.nodes()
        for node in all_nodes:
            if node == 'Input':
                colors.append('lightgreen')
                continue
            state = self.state[node]['state']
            if state == 'success':
                colors.append('lightgreen')
            elif state == 'failed':
                colors.append('orange')
            elif state == 'waiting':
                colors.append('lightgray')
            else:
                colors.append('lightblue')
        return dict(zip(all_nodes, colors))

    def get_label_dict(self):
        node_label_dict = dict()
        for each in self.graph.nodes():
            if each == 'Input':
                node_label_dict[each] = each
                continue
            used_time = self.state[each]['used_time']
            if isinstance(used_time, str):
                if used_time == 'unknown':
                    node_label_dict[each] = ''
                else:
                    node_label_dict[each] = used_time
            elif float(used_time) <= 0:
                node_label_dict[each] = ''
            else:
                node_label_dict[each] = str(used_time) + 's'
            node_label_dict[each] = each + '\n' + node_label_dict[each]
        return node_label_dict

    def draw(self, img_file='state.png'):
        self.add_edges()
        # pos = nx.kamada_kawai_layout(self.graph)
        # pos = nx.spring_layout(self.graph)
        pos = nx.nx_pydot.pydot_layout(self.graph, prog='dot')
        node_label_dict = self.get_label_dict()
        color_dict = self.get_color_dict()
        tmp_dict = dict()
        for k, v in color_dict.items():
            tmp_dict.setdefault(v, list())
            tmp_dict[v].append(k)
        plt.figure(figsize=(12, 8))
        for color, group in tmp_dict.items():
            if group[0] == 'Input':
                state = 'success'
            else:
                state = self.state[group[0]]['state']
            nx.draw(self.graph, pos=pos, nodelist=group, labels=node_label_dict,
                    with_labels=True, font_size=9, node_shape='o', node_size=1000,
                    node_color=color, label=state, alpha=1, width=0.7, style='dotted')

        # pos_attrs = {}
        # for node, coords in pos.items():
        #     pos_attrs[node] = (coords[0], coords[1] - 0.1)
        #
        # nx.draw_networkx_labels(self.graph, pos=pos_attrs, labels=node_label_dict,
        #                         font_size=8, alpha=0.8)
        plt.axis('off')
        plt.legend(loc='best', fontsize='small', markerscale=0.7, frameon=False)
        plt.savefig(img_file, dpi=200, bbox_inches='tight')
        plt.close()