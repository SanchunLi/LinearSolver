#include"bac.h"

int main(int argc, char** argv)
{

	clock_t start, end, first_start, first_end;
	first_start = clock();

	//读入case
	if (argc != 3)
		return false;
	string input_file(argv[1]);
	string output_file(argv[2]);

	BranchAndCut bac(input_file, output_file);

	Matf mat_a_colume = bac.store_by_colume(bac._mat_a);
	Matf mat_b_colume = bac.store_by_colume(bac._mat_b);

	double sum_demand = 0;
	//获取总的流量需求
	for (auto i : mat_b_colume[2])
	{
		sum_demand += i;
	}

	Vectorf node_max_flow = bac.one_node_max_flow(mat_a_colume); //
	int min_servers_number = bac.min_server_number(node_max_flow, sum_demand); //
	//int min_servers_number = ;

	//初始化A,cj,b
	//约束的数量
	/*对所用的数据结构先进行初始化，为了减少后期进行push_back时，数据移动的时间消耗*/
	int constraint_number = 2 * (bac.edge_number_ + bac.net_node_number_) + 1; //m 

	Matf initial_A(constraint_number);
	Vectorf initial_b(constraint_number, 0);
	Vectorf initial_cj(1 + 3 * bac.net_node_number_ + 4 * bac.edge_number_, 0); //n+m

	// 用于优化
	Mati initial_A_nonzero_mat(constraint_number);
	Vectori initial_cj_index(1 + 3 * bac.net_node_number_ + 4 * bac.edge_number_, 1);
	Vectori initial_cb_index(constraint_number, 1);

	//构造A
	bac.construct_A(initial_cj, initial_A, initial_b, node_max_flow, min_servers_number, constraint_number, 
		mat_a_colume, mat_b_colume,
		initial_A_nonzero_mat);



	//初始化上下界
	//注意，大规模问题9999999可能不够

	pair<double, Vectorf> gub;
	gub.first = bac.consumer_node_number_*bac.price_server_;
	//gub.first = bac.net_node_number_*bac.price_server_;
	//gub.first = 22000;



	//初始化问题表,并添加第一个节点
	map<double, vector<double>, greater<double >> node_map;
	//map<double, vector<double>> node_map;
	//Vectorf first_x(bac.consumer_node_number_, -2);
	//node_map.insert({ 0, first_x });


	double time = 0;

	//建立cj A b
	Matf A1(initial_A);
	Vectorf cj1(initial_cj), b1(initial_b);
	Solve first_node, first_solve, first_solve_erase_edge;
	//first_node.second = first_x;

	// first_node.first存放目标值  first.second解向量（只包含节点）(可行解时会再加入边)，  -2代表该节点未确定
	for (int i = 0; i < bac.net_node_number_; i++)
		first_node.second.push_back(-2);

	//进行对偶单纯形标准化和计算
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

	//根据计算结果是否是可信可行解，如果不是可行解，将里面最大的一个小数值取-1，剩下的遍历全为-2;
	//这样是为了分支时方便，这里只考虑服务器位置，因为流量肯定是整数的
	if (!bac.is_feasible_solution(first_solve.second))
	{
		//在解中去掉边
		first_solve_erase_edge.first = first_solve.first;
		int solve_size = first_solve.second.size();
		first_solve_erase_edge.second.reserve(solve_size - 2 * bac.edge_number_);  //
		
		
		copy(first_solve.second.begin(), first_solve.second.begin() + (solve_size - 2 * bac.edge_number_), back_inserter(first_solve_erase_edge.second));

		//在表示服务器的变量中找到最大的值 (该变量都是小数值)
		auto iter1 = max_element(first_solve_erase_edge.second.begin(), first_solve_erase_edge.second.end());
		*iter1 = -1;  //代表下一个要分支的节点
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


		//若node_map为空则说明遍历完，最优解已经有了,否则继续
		if (node_map.empty())
		{
			cout << "the opitimum has computed!" << endl;
			break;
		}

		//selection best-selection
		//在node_map中找到最大的一个，也就是最大的llb的节点分布进行计算,然后删除这个node
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

			
			//// first_node.first存放目标值  first.second解向量（只包含节点）(可行解时会再加入边)，  -2代表该节点未确定
			///*while()循环外面也有dualsimplex和判断是否是可行解, 是为了进入while(), 因为node_map.empty()*/ 
			//for (int i = 0; i < bac.net_node_number_; i++)
			//	first_node.second.push_back(-2);

			//进行对偶单纯形标准化和计算
			//bac.simplex_standard(cj, A, b);
			Mati A_nonzero_index(initial_A_nonzero_mat);
			Vectori cj_index(initial_cj_index), cb_index(initial_cb_index);

			//bac.dualsimplex(cj, A, b, first_solve, 0, A_nonzero_index_temp,cb_index,cj_index,0);

			//建立cj A b
			Matf A = initial_A;
			Vectorf cj = initial_cj;
			Vectorf b = initial_b;
			Solve solve;

			//根据选定的node，确定的值固定，也就是修改CJ,A，B中的行列值，这里的修改fix函数还有问题
			//再进行对偶单纯形的计算

			start = clock();
			//bac.fix(select_0_node.second, cj, A, b);
			int fixed_number = bac.fixed(cj, A, b, select_0_node.second, A_nonzero_index, cb_index, cj_index);
			//Mati fixed_Anonzero = bac.fix_A_nonzero_colume(select_0_node.second);
			//bac.simplex_standard(cj, A, b);
			bac.dualsimplex(cj, A, b, solve, 0, A_nonzero_index, cb_index, cj_index, fixed_number);

			++dual_count;
			end = clock();
			cout << end - start << endl;

			//非可行解退出循环
			if (solve.first == -1)
				continue;

			//根据计算结果是否是可信可行解，如果不是可行解，将里面最大的一个小数值取-1，剩下的遍历全为-2;
			//这样是为了分支时方便，这里只考虑服务器位置，因为流量肯定是整数的
			if (!bac.is_feasible_solution(solve.second))
			{

				//如果是非可行解,则进行对所求解中，最大小数标记，然后对这位数子进行分支成0 或 1.得到两个新的解，放入node_map数据结构中<double,vector<double>>

				//在解中去掉去掉边
				solve_erase_edge.first = solve.first;
				int solve_size = solve.second.size();
				solve_erase_edge.second.reserve(solve_size - 2 * bac.edge_number_);
				copy(solve.second.begin(), solve.second.begin() + (solve_size - 2 * bac.edge_number_), back_inserter(solve_erase_edge.second));

				//选取最大的的小数进行标记-1，为了下次branch

				auto iter1 = max_element(solve_erase_edge.second.begin(), solve_erase_edge.second.end());
				*iter1 = -1;
				//其余位置标记为-2，意思为自由变量
				for (int i = 0; i < solve_erase_edge.second.size(); i++)
				{
					if (solve_erase_edge.second[i] != -1)
						solve_erase_edge.second[i] = -2;
				}

				//根据服务器位置的解，计算出最终花费，并且将下一个分支点标为-1.存入map_node
				bac.recovery_x(solve_erase_edge, select_0_node);
				if (select_0_node.first < gub.first)
				{
					node_map.insert(select_0_node);
				}
			}/*end if not feasible solution*/
			else
			{
				//如果算出的结果是带边变量的，则先切除边，在计算最小花费。然后和全局上界进行对比。

				if (solve.second.size() > bac.edge_number_)
				{
					solve_erase_edge.first = solve.first;
					int solve_size = solve.second.size();
					solve_erase_edge.second.reserve(solve_size - 2 * bac.edge_number_);
					copy(solve.second.begin(), solve.second.begin() + (solve_size - 2 * bac.edge_number_), back_inserter(solve_erase_edge.second));

					bac.recovery_x(solve_erase_edge, select_0_node);

					//solve中还有服务器变量
					//加入边
					copy(solve.second.begin() + solve_size - 2 * bac.edge_number_, solve.second.end(), back_inserter(select_0_node.second));
					if (select_0_node.first < gub.first)
						gub = select_0_node;
				}
				else
				{
					//这个else是针对第一次solve没边的情况。直接计算花费，并和gub进行比较

					//************************************************************
					// 等到顶点全都(fixed)确定时（放或不放服务器）， solve里的解向量只有边变量，此时就不需要recovery_x函数。
					//只要固定了所有顶点，都会执行这一块代码
					select_0_node.first = 0;
					for (auto &i : select_0_node.second)
					{
						if (i == 1)
							select_0_node.first += bac.price_server_;
					}
					select_0_node.first += -solve.first;  //(solve.first只计算了边的流量费用),此时加上服务器费用
					copy(solve.second.begin(), solve.second.end(), back_inserter(select_0_node.second));//select_0_node只有顶点
					//solve里此时只有边变量
					if (select_0_node.first < gub.first)
						gub = select_0_node; //解向量中只有服务器变量（没有边变量）

				}

			}/*end else if(feasible solution)*/

			//	A.clear(),clear()只是析构，但是不释放内存
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
