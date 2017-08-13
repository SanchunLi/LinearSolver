#include"bac.h"

int main(int argc, char** argv)
{

	clock_t start, end, first_start, first_end;
	first_start = clock();

	//����case
	if (argc != 3)
		return false;
	string input_file(argv[1]);
	string output_file(argv[2]);

	BranchAndCut bac(input_file, output_file);

	Matf mat_a_colume = bac.store_by_colume(bac._mat_a);
	Matf mat_b_colume = bac.store_by_colume(bac._mat_b);

	double sum_demand = 0;
	//��ȡ�ܵ���������
	for (auto i : mat_b_colume[2])
	{
		sum_demand += i;
	}

	Vectorf node_max_flow = bac.one_node_max_flow(mat_a_colume); //
	int min_servers_number = bac.min_server_number(node_max_flow, sum_demand); //
	//int min_servers_number = ;

	//��ʼ��A,cj,b
	//Լ��������
	/*�����õ����ݽṹ�Ƚ��г�ʼ����Ϊ�˼��ٺ��ڽ���push_backʱ�������ƶ���ʱ������*/
	int constraint_number = 2 * (bac.edge_number_ + bac.net_node_number_) + 1; //m 

	Matf initial_A(constraint_number);
	Vectorf initial_b(constraint_number, 0);
	Vectorf initial_cj(1 + 3 * bac.net_node_number_ + 4 * bac.edge_number_, 0); //n+m

	// �����Ż�
	Mati initial_A_nonzero_mat(constraint_number);
	Vectori initial_cj_index(1 + 3 * bac.net_node_number_ + 4 * bac.edge_number_, 1);
	Vectori initial_cb_index(constraint_number, 1);

	//����A
	bac.construct_A(initial_cj, initial_A, initial_b, node_max_flow, min_servers_number, constraint_number, 
		mat_a_colume, mat_b_colume,
		initial_A_nonzero_mat);



	//��ʼ�����½�
	//ע�⣬���ģ����9999999���ܲ���

	pair<double, Vectorf> gub;
	gub.first = bac.consumer_node_number_*bac.price_server_;
	//gub.first = bac.net_node_number_*bac.price_server_;
	//gub.first = 22000;



	//��ʼ�������,����ӵ�һ���ڵ�
	map<double, vector<double>, greater<double >> node_map;
	//map<double, vector<double>> node_map;
	//Vectorf first_x(bac.consumer_node_number_, -2);
	//node_map.insert({ 0, first_x });


	double time = 0;

	//����cj A b
	Matf A1(initial_A);
	Vectorf cj1(initial_cj), b1(initial_b);
	Solve first_node, first_solve, first_solve_erase_edge;
	//first_node.second = first_x;

	// first_node.first���Ŀ��ֵ  first.second��������ֻ�����ڵ㣩(���н�ʱ���ټ����)��  -2����ýڵ�δȷ��
	for (int i = 0; i < bac.net_node_number_; i++)
		first_node.second.push_back(-2);

	//���ж�ż�����α�׼���ͼ���
	//bac.simplex_standard(cj, A, b);
	Mati A_nonzero_index_temp(initial_A_nonzero_mat);
	Vectori cj_index1(initial_cj_index), cb_index1(initial_cb_index);

	Vectorf all_sure; //(bac.net_node_number_,0);

	bool net = 0;
	int count = 0;
	multimap<int, int, greater<int>> sort_demand;

	//	for (int i = 0; i < mat_b_colume[2].size(); i++){
	//		sort_demand.insert({ mat_b_colume[2][i],mat_b_colume[1][i] });
	//	}

	//	int service = bac.consumer_node_number_ * 1 / 2;
	//	
	//	
	//	for (int i = 0; i <service ; i++){
	//		auto iter = sort_demand.begin();
	//		all_sure[iter->second] = 1;
	//		sort_demand.erase(iter);
	//	}
	//	int fixed_first=bac.fixed(cj1, A1, b1, all_sure, A_nonzero_index_temp, cb_index1, cj_index1);

	bac.dualsimplex(cj1, A1, b1, first_solve, 0, A_nonzero_index_temp, cb_index1, cj_index1, 0);
	//end = clock();

	//���ݼ������Ƿ��ǿ��ſ��н⣬������ǿ��н⣬����������һ��С��ֵȡ-1��ʣ�µı���ȫΪ-2;
	//������Ϊ�˷�֧ʱ���㣬����ֻ���Ƿ�����λ�ã���Ϊ�����϶���������
	if (!bac.is_feasible_solution(first_solve.second))
	{
		//�ڽ���ȥ����
		first_solve_erase_edge.first = first_solve.first;
		int solve_size = first_solve.second.size();
		first_solve_erase_edge.second.reserve(solve_size - 2 * bac.edge_number_);  //
		
		
		copy(first_solve.second.begin(), first_solve.second.begin() + (solve_size - 2 * bac.edge_number_), back_inserter(first_solve_erase_edge.second));

		//�ڱ�ʾ�������ı������ҵ�����ֵ (�ñ�������С��ֵ)
		auto iter1 = max_element(first_solve_erase_edge.second.begin(), first_solve_erase_edge.second.end());
		*iter1 = -1;  //������һ��Ҫ��֧�Ľڵ�
		for (int i = 0; i < bac.net_node_number_; i++)
		{
			if (first_solve_erase_edge.second[i] != -1)
				first_solve_erase_edge.second[i] = -2;
		}
		bac.recovery_x(first_solve_erase_edge, first_node);
		node_map.insert(first_node);
	}
	else
	{
		if (first_solve.first < gub.first)
			gub = first_solve;

	}

	//Matf(A1).swap(A1);
	//Vectorf(cj1).swap(cj1);
	//Vectorf(b1).swap(b1);
	Matf().swap(A1);
	Vectorf().swap(cj1);
	Vectorf().swap(b1);
	Mati().swap(A_nonzero_index_temp);
	Vectori().swap(cj_index1);
	Vectori().swap(cb_index1);

	int dual_count = 0;

	//***************************************************************************************
	//branch and cut
	// time is for time-limited
	while (!node_map.empty() && time<188000)
	{


		//��node_mapΪ����˵�������꣬���Ž��Ѿ�����,�������
		if (node_map.empty())
		{
			cout << "the opitimum has computed!" << endl;
			break;
		}

		//selection best-selection
		//��node_map���ҵ�����һ����Ҳ��������llb�Ľڵ�ֲ����м���,Ȼ��ɾ�����node
		auto iter = node_map.begin();
		Solve select_node = *iter;
		cout << "relaxtion feasible: " << select_node.first << endl;
		node_map.erase(iter);

		//branch
		for (int is_server = 0; is_server < 2; is_server++)
		{

			Solve select_0_node = select_node, solve_erase_edge;

			if (select_0_node.first < gub.first)
			{
				auto iter = find(select_0_node.second.begin(), select_0_node.second.end(), -1);
				*iter = is_server;
			}
			else
				continue;

			
			//// first_node.first���Ŀ��ֵ  first.second��������ֻ�����ڵ㣩(���н�ʱ���ټ����)��  -2����ýڵ�δȷ��
			///*while()ѭ������Ҳ��dualsimplex���ж��Ƿ��ǿ��н�, ��Ϊ�˽���while(), ��Ϊnode_map.empty()*/ 
			//for (int i = 0; i < bac.net_node_number_; i++)
			//	first_node.second.push_back(-2);

			//���ж�ż�����α�׼���ͼ���
			//bac.simplex_standard(cj, A, b);
			Mati A_nonzero_index(initial_A_nonzero_mat);
			Vectori cj_index(initial_cj_index), cb_index(initial_cb_index);

			//bac.dualsimplex(cj, A, b, first_solve, 0, A_nonzero_index_temp,cb_index,cj_index,0);

			//����cj A b
			Matf A = initial_A;
			Vectorf cj = initial_cj;
			Vectorf b = initial_b;
			Solve solve;

			//����ѡ����node��ȷ����ֵ�̶���Ҳ�����޸�CJ,A��B�е�����ֵ��������޸�fix������������
			//�ٽ��ж�ż�����εļ���

			start = clock();
			//bac.fix(select_0_node.second, cj, A, b);
			int fixed_number = bac.fixed(cj, A, b, select_0_node.second, A_nonzero_index, cb_index, cj_index);
			//Mati fixed_Anonzero = bac.fix_A_nonzero_colume(select_0_node.second);
			//bac.simplex_standard(cj, A, b);
			bac.dualsimplex(cj, A, b, solve, 0, A_nonzero_index, cb_index, cj_index, fixed_number);

			++dual_count;
			end = clock();
			cout << end - start << endl;

			//�ǿ��н��˳�ѭ��
			if (solve.first == -1)
				continue;

			//���ݼ������Ƿ��ǿ��ſ��н⣬������ǿ��н⣬����������һ��С��ֵȡ-1��ʣ�µı���ȫΪ-2;
			//������Ϊ�˷�֧ʱ���㣬����ֻ���Ƿ�����λ�ã���Ϊ�����϶���������
			if (!bac.is_feasible_solution(solve.second))
			{

				//����Ƿǿ��н�,����ж�������У����С����ǣ�Ȼ�����λ���ӽ��з�֧��0 �� 1.�õ������µĽ⣬����node_map���ݽṹ��<double,vector<double>>

				//�ڽ���ȥ��ȥ����
				solve_erase_edge.first = solve.first;
				int solve_size = solve.second.size();
				solve_erase_edge.second.reserve(solve_size - 2 * bac.edge_number_);
				copy(solve.second.begin(), solve.second.begin() + (solve_size - 2 * bac.edge_number_), back_inserter(solve_erase_edge.second));

				//ѡȡ���ĵ�С�����б��-1��Ϊ���´�branch

				auto iter1 = max_element(solve_erase_edge.second.begin(), solve_erase_edge.second.end());
				*iter1 = -1;
				//����λ�ñ��Ϊ-2����˼Ϊ���ɱ���
				for (int i = 0; i < solve_erase_edge.second.size(); i++)
				{
					if (solve_erase_edge.second[i] != -1)
						solve_erase_edge.second[i] = -2;
				}

				//���ݷ�����λ�õĽ⣬��������ջ��ѣ����ҽ���һ����֧���Ϊ-1.����map_node
				bac.recovery_x(solve_erase_edge, select_0_node);
				if (select_0_node.first < gub.first)
				{
					node_map.insert(select_0_node);
				}
			}/*end if not feasible solution*/
			else
			{
				//�������Ľ���Ǵ��߱����ģ������г��ߣ��ڼ�����С���ѡ�Ȼ���ȫ���Ͻ���жԱȡ�

				if (solve.second.size() > bac.edge_number_)
				{
					solve_erase_edge.first = solve.first;
					int solve_size = solve.second.size();
					solve_erase_edge.second.reserve(solve_size - 2 * bac.edge_number_);
					copy(solve.second.begin(), solve.second.begin() + (solve_size - 2 * bac.edge_number_), back_inserter(solve_erase_edge.second));

					bac.recovery_x(solve_erase_edge, select_0_node);

					//solve�л��з���������
					//�����
					copy(solve.second.begin() + solve_size - 2 * bac.edge_number_, solve.second.end(), back_inserter(select_0_node.second));
					if (select_0_node.first < gub.first)
						gub = select_0_node;
				}
				else
				{
					//���else����Ե�һ��solveû�ߵ������ֱ�Ӽ��㻨�ѣ�����gub���бȽ�

					//************************************************************
					// �ȵ�����ȫ��(fixed)ȷ��ʱ���Ż򲻷ŷ��������� solve��Ľ�����ֻ�б߱�������ʱ�Ͳ���Ҫrecovery_x������
					//ֻҪ�̶������ж��㣬����ִ����һ�����
					select_0_node.first = 0;
					for (auto &i : select_0_node.second)
					{
						if (i == 1)
							select_0_node.first += bac.price_server_;
					}
					select_0_node.first += -solve.first;  //(solve.firstֻ�����˱ߵ���������),��ʱ���Ϸ���������
					copy(solve.second.begin(), solve.second.end(), back_inserter(select_0_node.second));//select_0_nodeֻ�ж���
					//solve���ʱֻ�б߱���
					if (select_0_node.first < gub.first)
						gub = select_0_node; //��������ֻ�з�����������û�б߱�����

				}

			}/*end else if(feasible solution)*/

			//	A.clear(),clear()ֻ�����������ǲ��ͷ��ڴ�
			Matf().swap(A);
			Vectorf().swap(cj);
			Vectorf().swap(b);
			Mati().swap(A_nonzero_index);
			Vectori().swap(cj_index);
			Vectori().swap(cb_index);

		}/*end for() for branch*/

		first_end = clock();

		time = (first_end - first_start);

		
	}/*end while()*/

	if (!gub.second.empty())
	{
		copy(gub.second.begin(), gub.second.end(), back_inserter(all_sure));
		bac.write(all_sure, bac._output_file, mat_a_colume, mat_b_colume);
	}

	cout << "The best price is " << gub.first << endl;
	cout << "The dual number:  " << dual_count << endl;
	cout << time << "ms" << endl;
	system("pause");
	return 1;

}
