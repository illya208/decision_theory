#include <iostream>
#include <string>
#include <string.h>
#include <chrono>
#include<sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <algorithm>

#include <algorithm>
#include <functional>
#include <vector>
#include "windows.h"
#include "psapi.h"
#include <iomanip>
#define N 5
using namespace std;

void display(vector<int> &v)
{
	for (int i = 0; i < v.size(); i++)
	{
		cout << v[i] << " ";
	}
	cout << "\n" << endl;
}
class relation
{
public:
	relation() {

	};

	relation(const relation &copied) {
	}
	relation & operator=(const relation & copied) {}

	virtual void relation_union(relation rel) {};
	virtual void relation_intersection() {};
	virtual void relation_addition() {};
	virtual void relation_composition() {};
	virtual void relation_difference() {};
	virtual void relation_double_natured() {};
	virtual void relation_inverted() {};
	virtual void relation_narrowing() {};
	virtual void relation_symmetric_difference() {};
	virtual void relation_included() {};
	virtual void show() {};
	~relation() {};
private:

};


class relation_array : public virtual relation
{
public:
	relation_array() {
		Array = new  int*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new int[n];
	};

	relation_array(int m1, int n1) {
		m = m1;
		n = n1;
		Array = new  int*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new int[n];
	}

	relation_array(int m1, int n1, int** arr) {
		m = m1;
		n = n1;
		Array = new  int*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new int[n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				Array[i][j] = 0;
				Array[i][j] = arr[i][j];
			}
		}
	}

	relation_array(int m1, int n1, string elementary_level) { //empty,diagonal,full
		m = m1;
		n = n1;
		Array = new  int*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new int[n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (elementary_level.compare("input") == 0) {
					cin >> Array[i][j];
				}
				if (elementary_level.compare("empty") == 0) {
					Array[i][j] = 0;
				}
				if (elementary_level.compare("full") == 0) {
					Array[i][j] = 1;
				}
				if (elementary_level.compare("diagonal") == 0) {
					if (i == j) {
						Array[i][j] = 1;
					}
					else {
						Array[i][j] = 0;
					}
				}
				if (elementary_level.compare("antidiagonal") == 0) {
					if (i == j) {
						Array[i][j] = 0;
					}
					else {
						Array[i][j] = 1;
					}
				}
			}
		}
	}

	void show() {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout << Array[i][j] << " ";
			}cout << endl;
		}
		cout << endl;
	}

	relation_array(const relation_array & copied) { //copy contructor
		m = copied.m;
		n = copied.n;
		Array = new int*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new int[n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				Array[i][j] = copied.Array[i][j];
			}
		}
	}

	void relation_intersection(relation_array arr2) {
		relation_array result_arr(this->m, this->n);
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_intersection" << endl;
		}
		else {
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					if (this->Element(i, j) == 1 && arr2.Element(i, j) == 1) {
						result_arr.Array[i][j] = 1;
					}
					else {
						result_arr.Array[i][j] = 0;
					}
				}
			}
			result_arr.show();
		}
		*this = result_arr;
	}

	void relation_union(relation_array arr2) {
		relation_array result_arr(this->m, this->n);
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_union" << endl;
		}
		else {
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					if (this->Element(i, j) == 1 || arr2.Element(i, j) == 1) {
						result_arr.Array[i][j] = 1;
					}
					else {
						result_arr.Array[i][j] = 0;
					}
				}
			}
		}
		*this = result_arr;
	}

	void relation_difference(relation_array arr2) {
		relation_array result_arr(this->m, this->n);
		this->show();
		arr2.show();
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_difference" << endl;
		}
		else {
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					if (this->Element(i, j) == 1 && arr2.Element(i, j) == 0) {
						result_arr.Array[i][j] = 1;
					}
					else {
						result_arr.Array[i][j] = 0;
					}
				}
			}
		}
		*this = result_arr;

	}

	void relation_symmetric_difference(relation_array arr2) {
		relation_array result_arr(this->m, this->n);
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_difference" << endl;
		}
		else {
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					if (this->Element(i, j) == 1 && arr2.Element(i, j) == 0 || this->Element(i, j) == 0 && arr2.Element(i, j) == 1) {
						result_arr.Array[i][j] = 1;
					}
					else {
						result_arr.Array[i][j] = 0;
					}
				}
			}
		}
		*this = result_arr;
	}

	void relation_addition() {
		relation_array result_arr(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (this->Element(i, j) == 1) {
					result_arr.Array[i][j] = 0;
				}
				else {
					result_arr.Array[i][j] = 1;
				}
			}
		}
		*this = result_arr;
	}

	void relation_inverted() {
		relation_array result_arr(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				result_arr.Array[i][j] = this->Array[j][i];
			}
		}
		*this = result_arr;
	}

	void relation_composition(relation_array arr2) {
		relation_array result_arr(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				result_arr.Array[i][j] = 0;
			}
		}
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (this->Array[i][j] == 1) {
					for (int k = 0; k < this->n; k++) {
						if (arr2.Array[j][k] == 1) {
							result_arr.Array[i][k] = 1;
						}
					}
				}
			}
		}
		*this = result_arr;
	}

	void relation_narrowing(int position_x_1, int position_x_2) {
		int temp_m = m - 2;
		int temp_n = n - 2;
		relation_array result_arr(temp_m, temp_n, "empty");
		int tmpi = 0, tmpj = 0;
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (i != position_x_1 && i != position_x_2)
				{
					if (j != position_x_1 && j != position_x_2)
					{

						//cout<< tmpi<<" "<< tmpj<<" " << i<<" "<<j<<" "<< position_x_1 <<"  "<< position_x_2<<" " << this->Array[i][j] << endl;
						result_arr.Array[tmpi][tmpj] = this->Array[i][j];
						//cout << tmpi << " " << tmpj << " " << this->Array[i][j]<< endl;

						tmpj++;
						if (tmpj == n - 2) {
							tmpj = 0;
							tmpi++;
						}
					}
				}
			}

		}
		*this = result_arr;
	}

	void relation_narrowing(vector <int> position) {
		if (position.empty() != true) {
			int temp_m = m - 2;
			int temp_n = n - 2;
			relation_array result_arr(temp_m, temp_n, "empty");
			int tmpi = 0, tmpj = 0;
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					if (std::find(position.begin(), position.end(), i) == position.end()) {
						if (std::find(position.begin(), position.end(), j) == position.end()) {
							//cout << tmpi << " " << tmpj << " " << i << " " << j << " " << " " << this->Array[i][j] << endl;
							result_arr.Array[tmpi][tmpj] = this->Array[i][j]; // error
							//cout << tmpi << " " << tmpj << " " << this->Array[i][j]<< endl;

							tmpj++;
							if (tmpj == n - 2) {
								tmpj = 0;
								tmpi++;
							}
						}
					}

				}

			}
			*this = result_arr;
		}

	}




	void relation_double_natured() {
		relation_array result_arr(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (this->Element(i, j) == 1) {
					result_arr.Array[i][j] = 0;
				}
				else {
					result_arr.Array[i][j] = 1;
				}
			}
		}
		result_arr.relation_inverted();
		*this = result_arr;
	}

	bool relation_included(relation_array arr2) {
		bool included_array = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (arr2.Array[i][j] == 1 && Array[i][j] == 0) {
					included_array = false;
					break;
				}
				/*else {
					included_array = true;
				}*/
			}
		}
		/*if (included_array == true) {
			cout << "included" << endl;
		}
		else {
			cout << "not included" << endl;
		}*/
		return included_array;
	}

	void calculate(relation_array Q, relation_array R) { //this==P		
		auto  now = std::chrono::high_resolution_clock::now();
		this->relation_composition(Q);
		R.relation_double_natured();
		this->relation_difference(R);
		auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
		std::cout << "Time : " << elapsed.count() << "ns.\n";

	}
	///
	///
	///
	//2 laba
	bool check_Reflexive() {
		bool reflexibility = true;
		for (int i = 0; i < m; i++) {
			if (this->Array[i][i] == 1) {
				reflexibility = true;
			}
			else {
				reflexibility = false;
				break;
			}
		}
		if (reflexibility == true) {
			cout << "relation Reflexive" << endl;
		}
		else {
			cout << "relation is not Reflexive" << endl;
		}
		return reflexibility;

	}

	void check_antireflexive() {
		bool antireflexibility = true;
		for (int i = 0; i < m; i++) {
			if (this->Array[i][i] == 0) {
				antireflexibility = true;
			}
			else {
				antireflexibility = false;
				break;
			}
		}
		if (antireflexibility) {
			cout << "relation antireflexive" << endl;
		}
		else {
			cout << "relation is not antireflexive" << endl;
		}
	}
	bool check_Symmetric() {
		bool symetric = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array[i][j] == Array[j][i]) {
					symetric = true;
				}
				else {
					symetric = false;
					break;
				}
			}
		}
		if (symetric) {
			cout << "relation symetric" << endl;
		}
		else {
			cout << "relation is not symetric" << endl;
		}
		return symetric;
	}


	bool check_Asymmetric() {
		bool asymmetric = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array[i][i] == 0) {
					if (Array[i][j] == 0 && Array[j][i] == 1 || Array[i][j] == 1 && Array[j][i] == 0 || Array[i][j] == 0 && Array[j][i] == 0) {
						asymmetric = true;
					}
					else {
						asymmetric = false;
						break;
					}
				}
				else {
					asymmetric = false;
					break;
				}
			}
		}
		if (asymmetric) {
			cout << "relation asymmetric" << endl;
		}
		else {
			cout << "relation is not asymmetric" << endl;
		}
		return asymmetric;
	}


	bool check_Antisymmetric() {
		bool antisymmetric = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					if (Array[i][j] == 0 && Array[j][i] == 1 || Array[i][j] == 1 && Array[j][i] == 0 || Array[i][j] == 0 && Array[j][i] == 0) {
						antisymmetric = true;
					}
					else {
						antisymmetric = false;
						break;
					}
				}
			}
		}
		if (antisymmetric) {
			cout << "relation antisymmetric" << endl;
		}
		else {
			cout << "relation is not antisymmetric" << endl;
		}
		return antisymmetric;
	}

	bool check_transitivity() {
		relation_array copy(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copy.Array[i][j] = Array[i][j];
			}
		}

		copy.relation_composition(*this);

		bool inc = this->relation_included(copy);
		if (inc == true) {
			cout << "transitive" << endl;
		}
		else {
			cout << "not transitive" << endl;
		}
		return inc;
	}

	void check_Acyclic() {

		bool Acyclic = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array[i][j] == 1 && Array[j][i] == 1) {
					Acyclic = false;
					break;
				}
			}
		}
		if (Acyclic == true) {
			cout << "Acyclic" << endl;
		}
		else {
			cout << "not Acyclic" << endl;
		}
	}

	bool check_linear() { // the same as Зв’язний
		bool linear = true;
		relation_array copy(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copy.Array[i][j] = Array[i][j];
			}
		}
		relation_array copyw(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copyw.Array[i][j] = Array[i][j];
			}
		}

		copy.relation_inverted();
		this->relation_union(copy);
		relation_array E(this->m, this->n, "diagonal");
		relation_array U(this->m, this->n, "full");
		this->relation_difference(E);
		U.relation_difference(E);

		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (this->Array[i][j] != U.Array[i][j]) {
					linear = false;
					break;
				}
			}
		}
		*this = copyw;
		return linear;
	}


	//2толерантності, еквівалентності, квазіпорядку, порядку, строгого порядку, лінійного
	//порядку, строгого лінійного порядку.

	bool check_tolerance() {
		bool tolerance;
		if (this->check_Reflexive() == true && this->check_Symmetric() == true)
		{
			tolerance = true;
			cout << "tolerance" << endl;
		}
		else {
			tolerance = false;
		}
		return tolerance;
	}

	bool check_equivalence() {
		bool equivalence;
		if (this->check_Reflexive() == true && this->check_Symmetric() == true && this->check_transitivity() == true)
		{
			equivalence = true;
			cout << "equivalence" << endl;
		}
		else {
			equivalence = false;
		}
		return equivalence;
	}

	bool check_quasi_order() {
		bool quasi;
		if (this->check_Reflexive() == true && this->check_transitivity() == true)
		{
			quasi = true;

			cout << "quasi_order" << endl;
		}
		else {
			quasi = false;
		}
		return quasi;
	}

	bool check_order() {
		bool order;
		if (this->check_Reflexive() == true && this->check_Antisymmetric() == true && this->check_transitivity() == true)
		{
			order = true;

			cout << "order" << endl;
		}
		else {
			order = false;
		}
		return order;
	}

	bool check_strict_order() {
		bool strict_order;
		if (this->check_Asymmetric() == true && this->check_transitivity() == true)
		{
			strict_order = true;

			cout << "strict_order" << endl;
		}
		else {
			strict_order = false;
		}
		return strict_order;
	}
	bool check_linear_order() {
		bool linear_order;
		if (this->check_Reflexive() == true && this->check_Antisymmetric() == true
			&& this->check_transitivity() == true && this->check_linear() == true)
		{
			linear_order = true;

			cout << "linear_order" << endl;
		}
		else {
			linear_order = false;
		}
		return linear_order;
	}
	//TODO not known - not present in book
	bool check_strict_linear_order() {
		bool strict_linear_order;
		if (this->check_Asymmetric() == true && this->check_transitivity() == true && this->check_Antisymmetric()
			&& this->check_Reflexive() == true && this->check_transitivity() == true)
		{
			strict_linear_order = true;

			cout << "strict_linear_order" << endl;
		}
		else {
			strict_linear_order = false;
		}
		return strict_linear_order;

	}
	//help func for |
	bool is_equal(relation_array arr) {
		bool is_equal = true;
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (arr.Array[i][j] != Array[i][j]) {
					is_equal = false;
					break;
				}
			}
		}
		return is_equal;
	}

	bool check_vector(std::vector<relation_array> g1) {
		for (auto i = g1.begin(); i != g1.end(); ++i) {
			for (auto j = g1.begin(); j != g1.end(); ++j) {
				if (i->is_equal(*j) == true) {
					return true;
				}
			}
		}
	}
	//end of help func 




	//
	/*Для базового класу або його реалізацій написати функції знаходження симетричної та
		асиметричної складової, транзитивного замикання, досягальності та взаємної
		досягальності.*/
		//

	relation_array find_symmetric_component() {
		relation_array copy(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copy.Array[i][j] = Array[i][j];
			}
		}
		relation_array copy2(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copy2.Array[i][j] = Array[i][j];
			}
		}
		copy2.relation_inverted();
		copy.relation_intersection(copy2);
		return copy;

	}
	relation_array find_asymmetric_component() {
		relation_array copy(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copy.Array[i][j] = Array[i][j];
			}
		}
		relation_array copy2(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				copy2.Array[i][j] = Array[i][j];
			}
		}
		copy = copy.find_symmetric_component();
		copy2.relation_difference(copy);
		return copy2;
	}
	void find_transient_closure() {
		std::vector<relation_array> v;
		relation_array copy(*this);
		this->relation_composition(copy);//this is updated;
		v.push_back(*this);

		while (check_vector(v) != true) {
			this->relation_composition(copy);
			v.push_back(*this);
		}
		for (auto i = v.begin(); i != v.end(); ++i) {
			copy.relation_union(*i);
		}
		copy.show();
		*this = copy;

	}
	relation_array find_reach() {
		std::vector<relation_array> v;
		relation_array copy(*this);

		relation_array copy2(*this);
		this->relation_composition(copy);//this is updated;
		v.push_back(*this);

		while (check_vector(v) != true) {
			this->relation_composition(copy);
			v.push_back(*this);
		}
		for (auto i = v.begin(); i != v.end(); ++i) {
			copy.relation_union(*i);
		}
		copy.show();
		relation_array P(this->m, this->n, "diagonal");
		P.relation_union(copy);
		P.show();
		*this = copy2;
		return P;

	}
	relation_array find_mutual_reach() {

		relation_array copy(find_reach());
		relation_array copy2(copy);
		copy.relation_inverted();
		copy2.relation_intersection(copy);
		copy2.show();
		return copy2;
	}

	void check_type() {
		this->check_tolerance();
		cout << "!" << endl;
		this->check_equivalence();
		cout << "!" << endl;
		this->check_quasi_order();
		cout << "!" << endl;
		this->check_order();
		cout << "!" << endl;
		this->check_strict_order();
		cout << "!" << endl;
		this->check_linear_order();
		cout << "!" << endl;
		this->check_strict_linear_order();
		bool tr = this->check_transitivity();
		if (tr == false) {
			cout << "lets find transitive closure" << endl;
			this->find_transient_closure();
			this->check_tolerance();
			this->check_equivalence();
			this->check_quasi_order();
			this->check_order();
			this->check_strict_order();
			this->check_linear_order();
			this->check_strict_linear_order();
			bool tr = this->check_transitivity();
		}
	}









	//Для реалізацій базового класу написати функцію факторизації за заданим відношенням
	//еквівалентності.

	relation_array factorization_by_the_given_relation_equivalence(relation_array D) {
		std::vector<int> v;


		int ** a = new  int*[D.m];
		for (int z = 0; z < D.m; z++)
			a[z] = new int[D.n];

		for (int i = 0; i < D.m; i++) {
			for (int j = 0; j < D.n; j++) {
				a[i][j] = D.Array[i][j];
			}
		}
		int	m = D.m;
		int n = D.n;
		int rowTo = 0;
		size_t k = 0;

		for (size_t i = 0; i < m; i++)
		{
			size_t j = 0;

			while (j < k && not std::equal(a[i], a[i] + n, a[j])) j++;

			if (j == k)
			{
				if (k != i)
				{
					std::copy(a[i], a[i] + n, a[k]);
				}
				++k;
			}
		}

		if (k != m)
		{
			int **tmp = new int *[k];

			std::copy(a, a + k, tmp);
			std::for_each(a + k, a + m, std::default_delete<int[]>());

			delete[] a;

			a = tmp;
		}


		for (int i = 0; i < k; i++)
		{
			for (int j = 0; j < n; ++j) {
				cout << a[i][j] << " ";
			}cout << endl;
		}




		int ** m2 = new  int*[k];
		for (int z = 0; z < k; z++)
			m2[z] = new int[k];
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				m2[i][j] = 0;
			}
		}

		std::vector<std::vector<int> > imatrix;
		std::vector<int> w;
		for (int i = 0; i < k; i++) {
			for (int jj = 0; jj < n; jj++) {

				if (a[i][jj] == 1) {
					w.push_back(jj);
				}
			}
			for (auto ja = w.begin(); ja != w.end(); ja++) {
				cout << *ja << endl;
			}
			cout << "ss" << endl;
			imatrix.push_back(w);
			w.clear();

		}


		int point1 = 0, point2 = 0;
		for (vector<vector<int> >::iterator it = imatrix.begin(); it != imatrix.end(); it++) {
			for (vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); it2++) {
				for (vector<vector<int> >::iterator it3 = imatrix.begin(); it3 != imatrix.end(); it3++) {
					for (vector<int>::iterator it4 = (*it3).begin(); it4 != (*it3).end(); it4++) {
						cout << *it2 << "    " << *it4 << endl;
						if (this->Array[*it2][*it4] == 1) {
							m2[point1][point2] = 1;
							cout << *it2 << "    " << *it4 << "  point 1  " << point1 << "  point 2  " << point2 << endl;

							//point2++;
						}
						else {
							m2[point1][point2] = 0;
							//point2++;
						}

					}point2++;
				}/*oint2++;*/

			}
			cout << endl;
			point2 = 0;
			point1++;
		}
		cout << "K=====" << k << endl;
		for (int ic = 0; ic < k; ic++) {
			for (int jz = 0; jz < k; jz++) {
				cout << m2[ic][jz] << "  ";

			}cout << endl;
		}

		cout << "K=====" << k << endl;
		relation_array created(k, k, "empty");
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				created.Array[i][j] = m2[i][j];
			}
		}

		return created;


	}

	//test 4 laba
	void find_maximum() {
		bool flag = true;
		int count = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array[i][j] == 1) {
					flag = true;
				}
				else {
					flag = false;
					break;
				}

			}
			if (flag) {
				cout << "alternative " << i + 1 << " is maximum" << endl;
				count++;
			}
		}

		if (count == 0) {
			cout << "maximum not found" << endl;
		}
	}

	void find_majorant() {
		bool flag = true;
		int count = 0;

		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				if (Array[i][j] == 0) {
					flag = true;
				}
				else {
					flag = false;
					break;
				}

			}
			if (flag) {
				cout << "alternative " << j + 1 << " is majorant" << endl;
				count++;
			}
		}

		if (count == 0) {
			cout << "majorant not found" << endl;
		}
	}












	relation_array & operator=(const relation_array & copied) { // overloaded = 
		if (this != &copied) {
			//std::cout << "overloaded operator = "<< std::endl;
			for (int z = 0; z < m; z++) //clean buff
				delete[] Array[z];
			delete[] Array;
			m = copied.m;
			n = copied.n;
			Array = new int*[m];
			for (int z = 0; z < m; z++)
				Array[z] = new int[n];
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					Array[i][j] = copied.Array[i][j];
				}
			}
		}
		else {
			std::cout << "own=" << std::endl;
		}
		return *this;
	}

	int& Element(int i, int j)
	{
		return Array[i][j];
	}

	~relation_array() {
		for (int z = 0; z < m; z++)
			delete[] Array[z];
		delete[] Array;
	}

	int **Array;
	int m;
	int n;
};


class relation_cat : public virtual relation // 
{
public:
	relation_cat(int m1, int n1, string elementary_level) {
		this->m = m1;
		this->n = n1;
		string x = "x";
		Array = new  string*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new string[n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (elementary_level.compare("input") == 0) {
					cin >> Array[i][j];
				}
				if (elementary_level.compare("empty") == 0) {
					Array[i][j] = "0";
				}
				if (elementary_level.compare("full") == 0) {
					stringstream ss;
					int buf = j + 1;
					ss << buf;
					string str = ss.str();
					Array[i][j] = x + str;
				}
				if (elementary_level.compare("diagonal") == 0) {
					if (i == j) {
						stringstream ss;
						int buf = j + 1;
						ss << buf;
						string str = ss.str();
						Array[i][j] = x + str;
					}
					else {
						Array[i][j] = "0";
					}
				}
				if (elementary_level.compare("antidiagonal") == 0) {
					if (i == j) {
						Array[i][j] = "0";
					}
					else {
						stringstream ss;
						int buf = j + 1;
						ss << buf;
						string str = ss.str();
						Array[i][j] = x + str;
					}
				}
			}
		}
	}

	relation_cat(int m1, int n1) {
		m = m1;
		n = n1;
		Array = new  string*[m];
		for (int z = 0; z < m; z++)
			Array[z] = new string[n];
	}

	void show() {
		for (int i = 0; i < m; i++) {
			cout << "  { ";
			for (int j = 0; j < n; j++) {
				if (Array[i][j].compare("0") != 0) {
					cout << Array[i][j] << ", ";
				}
			}
			cout << " }  ";
		}
		cout << endl;
	}



	relation_cat from_matrix_to_cat(relation_array recovery) {
		relation_cat rc_res(recovery.m, recovery.n);
		string x = "x";
		recovery.relation_inverted();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (recovery.Array[i][j] == 1) {
					stringstream ss;
					int buf = j + 1;
					ss << buf;
					string str = ss.str();
					rc_res.Array[i][j] = x + str;
				}
				else {
					rc_res.Array[i][j] = "0";
				}
			}
		}
		return rc_res;
	}

	relation_array from_cat_to_matrix() {
		relation_array recovery(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array[i][j].compare("0") != 0) {
					recovery.Array[j][i] = 1;
				}
				else {
					recovery.Array[j][i] = 0;
				}
			}
		}
		return recovery;
	}

	void relation_union(relation_cat cat) {
		relation_array rc = this->from_cat_to_matrix();
		relation_array rc2 = cat.from_cat_to_matrix();
		rc.relation_union(rc2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_difference(relation_cat cat) {
		relation_array rc = this->from_cat_to_matrix();
		relation_array rc2 = cat.from_cat_to_matrix();
		rc.relation_difference(rc2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_intersection(relation_cat cat) {
		relation_array rc = this->from_cat_to_matrix();
		relation_array rc2 = cat.from_cat_to_matrix();
		rc.relation_intersection(rc2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_composition(relation_cat cat) {
		relation_array rc = this->from_cat_to_matrix();
		relation_array rc2 = cat.from_cat_to_matrix();
		rc.relation_composition(rc2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_double_natured() {
		relation_array rc = this->from_cat_to_matrix();
		rc.relation_double_natured();
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}


	void relation_inverted() {
		relation_array rc = this->from_cat_to_matrix();
		rc.relation_inverted();
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_addition() {
		relation_array rc = this->from_cat_to_matrix();
		rc.relation_addition();
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_narrowing(int position_x_1, int position_x_2) {
		relation_array rc = this->from_cat_to_matrix();
		rc.relation_narrowing(position_x_1, position_x_2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	void relation_symmetric_difference(relation_cat cat) {
		relation_array rc = this->from_cat_to_matrix();
		relation_array rc2 = cat.from_cat_to_matrix();
		rc.relation_symmetric_difference(rc2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}

	virtual void relation_included(relation_cat cat) {
		relation_array rc = this->from_cat_to_matrix();
		relation_array rc2 = cat.from_cat_to_matrix();
		rc.relation_included(rc2);
		relation_cat copy = from_matrix_to_cat(rc);
		*this = copy;
	}



	void calculate(relation_cat Q, relation_cat R) { //this==P		
		auto  now = std::chrono::high_resolution_clock::now();
		this->relation_composition(Q);
		R.relation_double_natured();
		this->relation_difference(R);
		auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
		std::cout << "Time : " << elapsed.count() << "ns.\n";

	}


	relation_cat & operator=(const relation_cat & copied) { // overloaded = 
		if (this != &copied) {
			for (int z = 0; z < m; z++) //clean buff
				delete[] Array[z];
			delete[] Array;
			m = copied.m;
			n = copied.n;
			Array = new string*[m];
			for (int z = 0; z < m; z++)
				Array[z] = new string[n];
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					Array[i][j] = copied.Array[i][j];
				}
			}
		}
		else {
			std::cout << "own=" << std::endl;
		}
		return *this;
	}

	string **Array;
	int m;
	int n;

private:
};





class relation_array_metr : public virtual relation_array {


public:
	/*relation_array_metr():relation_array {
		Array_metr = new  float*[m];
		for (int z = 0; z < m; z++)
			Array_metr[z] = new float[n];
	};*/

	relation_array_metr(int m1 = 5, int n1 = 5) :relation_array(m1, n1) {
		m = m1;
		n = n1;
		Array_metr = new  float*[m];
		for (int z = 0; z < m; z++)
			Array_metr[z] = new float[n];
	}

	relation_array_metr(int m1, int n1, float** arr) :relation_array(m1, n1) {
		m = m1;
		n = n1;
		Array_metr = new  float*[m];
		for (int z = 0; z < m; z++)
			Array_metr[z] = new float[n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				Array_metr[i][j] = arr[i][j];
			}
		}
	}
	relation_array_metr(int m1, int n1, string elementary_level) :relation_array(m1, n1) { //empty,full
		m = m1;
		n = n1;
		Array_metr = new  float*[m];
		for (int z = 0; z < m; z++)
			Array_metr[z] = new float[n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (elementary_level.compare("input") == 0) {
					cin >> Array_metr[i][j];
				}
				if (elementary_level.compare("empty") == 0) {
					Array_metr[i][j] = 0;
				}
				if (elementary_level.compare("full") == 0) {
					Array_metr[i][j] = 1;
				}
			}
		}
		/*for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if(Array_metr[i][j] !=0)
			}
		}*/
	}

	void relation_narrowing(vector <int> position) {
		if (position.empty() != true) {
			int temp_m = m - position.size();
			int temp_n = n - position.size();
			relation_array_metr result_arr(temp_m, temp_n, "empty");
			int tmpi = 0, tmpj = 0;


			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					for (int k = 0; k < position.size(); k++) {
						if (i == position[k] || j == position[k]) {
							break;
						}
						else {
							cout << this->Array_metr[i][j] << " ";
							result_arr.Array_metr[tmpi][tmpj] = this->Array_metr[i][j];
							tmpj++;
							if (tmpj == n - position.size()) {
								tmpj = 0;
								tmpi++;
							}

						}
					}
				}
				cout << endl;
			}











			//for (int i = 0; i < this->m; i++) {
			//	for (int j = 0; j < this->n; j++) {
			//		if (std::find(position.begin(), position.end(), i) == position.end()) {
			//			if (std::find(position.begin(), position.end(), j) == position.end()) {
			//				if (tmpj == n - position.size()) {
			//					tmpj = 0;
			//					tmpi++;
			//				}
			//				//cout << tmpi << " " << tmpj << " " << i << " " << j << " " << " " << this->Array[i][j] << endl;
			//				result_arr.Array_metr[tmpi][tmpj] = this->Array_metr[i][j]; // error
			//				//cout << tmpi << " " << tmpj << " " << this->Array[i][j]<< endl;

			//				tmpj++;
			//				
			//			}
			//		}

			//	}

			//}
			*this = result_arr;
		}

	}

	float get_v(int i, int j) {
		return Array_metr[i][j];
	}
	void set_v(int i, int j, float value) {
		Array_metr[i][j] = value;
	}


	bool help_func() {
		bool type = false;
		string type_str;
		bool is_once = false;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (Array_metr[i][j] == Array_metr[i][k] + Array_metr[k][j]) {
						//type = true;
						is_once = true;
						break;
					}
				}
				if (is_once == false) {
					type = false;
					return false;
					//type_str = "Multiplicative";
					break;
				}
				else {
					is_once = false;
					type = true;

				}
			}
		}
		return type;
	}

	string determining_type_of_relation() {
		bool type = this->help_func();
		string type_str;
		if (type == true) {
			type_str = "Additive";
		}
		else {
			type_str = "Multiplicative";
		}
		return type_str;
	}


	bool help_check() {

		bool match = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array_metr[i][j] == -1.0 * Array_metr[j][i] || (Array_metr[i][j] == 0 && Array_metr[j][i] == 0)) {
					//	cout << Array_metr[i][j] << " " << Array_metr[j][i] << " " << -1 * Array_metr[j][i] << endl;
					match = true;
				}
				else {
					match = false;
					return false;
					break;
				}

			}
		}
		return match;
	}

	string check_matched() {
		bool match = help_check();

		string type_str;
		if (match == true) {
			type_str = "Uzgodjena";
		}
		else {
			type_str = "Ne Uzgodjena";
		}
		return type_str;
	}


	void relation_union(relation_array_metr arr2) {
		relation_array_metr result_arr(this->m, this->n);
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_union" << endl;
		}
		else {

			if (this->determining_type_of_relation().compare("Additive") == 0 && arr2.determining_type_of_relation().compare("Additive") == 0) {
				for (int i = 0; i < this->m; i++) {
					for (int j = 0; j < this->n; j++) {
						if (Array_metr[i][j] == 0) {
							result_arr.Array_metr[i][j] = arr2.Array_metr[i][j];
						}
						else if (arr2.Array_metr[i][j] == 0) {
							result_arr.Array_metr[i][j] = Array_metr[i][j];
						}
						else {
							result_arr.Array_metr[i][j] = (Array_metr[i][j] + arr2.Array_metr[i][j]) / 2;
						}
					}
				}
			}
			else if (this->determining_type_of_relation().compare("Multiplicative") == 0 && arr2.determining_type_of_relation().compare("Multiplicative") == 0) {
				for (int i = 0; i < this->m; i++) {
					for (int j = 0; j < this->n; j++) {

						if (Array_metr[i][j] == 0) {
							result_arr.Array_metr[i][j] = arr2.Array_metr[i][j];
						}
						else if (arr2.Array_metr[i][j] == 0) {
							result_arr.Array_metr[i][j] = Array_metr[i][j];
						}
						else {
							result_arr.Array_metr[i][j] = pow((Array_metr[i][j] * arr2.Array_metr[i][j]), 0.5);
						}
					}
				}

			}
			else {
				cout << "oups - smth wrong" << endl;
			}

		}
		*this = result_arr;
	}

	void relation_intersection(relation_array_metr arr2) {
		relation_array_metr result_arr(this->m, this->n);
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_intersection" << endl;
		}
		else {

			if (this->determining_type_of_relation().compare("Additive") == 0 && arr2.determining_type_of_relation().compare("Additive") == 0) {
				for (int i = 0; i < this->m; i++) {
					for (int j = 0; j < this->n; j++) {
						if (Array_metr[i][j] != 0 && arr2.Array_metr[i][j] != 0) {
							result_arr.Array_metr[i][j] = (Array_metr[i][j] + arr2.Array_metr[i][j]) / 2;
						}
						else {
							result_arr.Array_metr[i][j] = 0;
						}
					}
				}
			}
			else if (this->determining_type_of_relation().compare("Multiplicative") == 0 && arr2.determining_type_of_relation().compare("Multiplicative") == 0) {
				for (int i = 0; i < this->m; i++) {
					for (int j = 0; j < this->n; j++) {
						if (Array_metr[i][j] != 0 && arr2.Array_metr[i][j] != 0) {
							result_arr.Array_metr[i][j] = pow((Array_metr[i][j] * arr2.Array_metr[i][j]), 0.5);
						}
						else {
							result_arr.Array_metr[i][j] = 0;
						}
					}
				}

			}

			else {
				cout << "oups - smth wrong" << endl;
			}

		}
		*this = result_arr;
	}

	void relation_difference(relation_array_metr arr2) {
		relation_array_metr result_arr(this->m, this->n);
		if (this->n != arr2.n || this->m != arr2.m) {
			cout << "error in relation_difference" << endl;
		}
		else {
			for (int i = 0; i < this->m; i++) {
				for (int j = 0; j < this->n; j++) {
					if (arr2.Array_metr[i][j] == 0) {
						result_arr.Array_metr[i][j] = Array_metr[i][j];
					}
					else {
						result_arr.Array_metr[i][j] = 0;
					}
				}
			}
		}
		*this = result_arr;

	}

	void check_transitive() {
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				if (this->Array_metr[i][j] != 0) {
					this->Array[i][j] = 1;
				}
				else {
					this->Array[i][j] = 0;
				}
			}
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout << Array[i][j] << "   ";
			}cout << endl;
		}
		cout << endl;



		this->check_transitivity();
	}

	void relation_composition(relation_array_metr arr2) {
		relation_array_metr result_arr(this->m, this->n);
		for (int i = 0; i < this->m; i++) {
			for (int j = 0; j < this->n; j++) {
				result_arr.Array_metr[i][j] = 0;
			}
		}
		result_arr.check_transitivity();
		if (this->determining_type_of_relation().compare("Additive") == 0 && arr2.determining_type_of_relation().compare("Additive") == 0) {
			cout << "Additive composition" << endl;
			std::vector<float> w;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
				{
					for (int v = 0; v < n; v++)
					{
						if (Array_metr[i][v] != 0 && arr2.Array_metr[v][j] != 0) {
							float num = this->Array_metr[i][v] + arr2.Array_metr[v][j];
							w.push_back(num);
						}
					}
					if (w.empty() == false) {
						float total = 0.0;
						for (auto i = w.begin(); i != w.end(); ++i) {
							total += *i;
						}
						total /= w.size();
						result_arr.Array_metr[i][j] = total;
					}
					else {
						result_arr.Array_metr[i][j] = 0;
					}
					w.clear();
				}
			}
		}
		else if (this->determining_type_of_relation().compare("Multiplicative") == 0 && arr2.determining_type_of_relation().compare("Multiplicative") == 0) {
			cout << "Multiplicative comosition" << endl;
			std::vector<float> w;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
				{
					for (int v = 0; v < n; v++)
					{
						if (Array_metr[i][v] != 0 && arr2.Array_metr[v][j] != 0) {
							float num = this->Array_metr[i][v] * arr2.Array_metr[v][j];
							w.push_back(num);
						}
					}
					if (w.empty() == false) {
						float total = 1.0;
						for (auto i = w.begin(); i != w.end(); ++i) {
							total *= *i;
						}
						total = pow(total, 1.0 / w.size());
						result_arr.Array_metr[i][j] = total;
					}
					else {
						result_arr.Array_metr[i][j] = 0;
					}
					w.clear();
				}
			}
		}
		else {
			cout << "oups - smth wrong" << endl;
		}

		*this = result_arr;
	}


	void show() {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout << left << fixed << setw(5) << setprecision(4) << Array_metr[i][j] << " ";
			}cout << endl;
		}
		cout << endl;
		cout << endl;
	}
	//4 laba



	void find_maximum() {
		bool flag = true;
		int count = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array_metr[i][j] == 1.0) {
					flag = true;
				}
				else {
					flag = false;
					break;
				}

			}
			if (flag) {
				cout << "alternative " << i + 1 << " is maximum" << endl;
				count++;
			}
		}

		if (count == 0) {
			cout << "maximum not found" << endl;
		}
	}

	void find_minimum() {
		bool flag = true;
		int count = 0;

		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				if (Array_metr[i][j] == 1) {
					flag = true;
				}
				else {
					flag = false;
					break;
				}

			}
			if (flag) {
				cout << "alternative " << j + 1 << " is minimum" << endl;
				count++;
			}
		}

		if (count == 0) {
			cout << "minimum not found" << endl;
		}

	}

	void find_majorant() {
		bool flag = true;
		int count = 0;

		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				if (Array_metr[i][j] == 0) {
					flag = true;
				}
				else {
					flag = false;
					break;
				}

			}
			if (flag) {
				cout << "alternative " << j + 1 << " is majorant" << endl;
				count++;
			}
		}

		if (count == 0) {
			cout << "majorant not found" << endl;
		}
	}

	void find_minorant() {
		bool flag = true;
		int count = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array_metr[i][j] == 0) {
					flag = true;
				}
				else {
					flag = false;
					break;
				}

			}
			if (flag) {
				cout << "alternative " << i + 1 << " is minorant" << endl;
				count++;
			}
		}

		if (count == 0) {
			cout << "minorant not found" << endl;
		}
	}




	void find_measure_of_proximity_equivalence(relation_array_metr Q2) {
		/*	vector<string> vector_of_Q1;
	*/
		vector<vector<int>> vector_of_Q1;
		vector<int> tmp1;
		for (int i = 0; i < m; i++) {
			bool flag = true;
			/*string names_of_row;*/
			for (int j = 0; j < n; j++) {
				if (Array_metr[i][j] == Array_metr[j][i] && Array_metr[i][j] == 1) {
					flag = true;
					/*	names_of_row += to_string(j);
							names_of_row += " ";*/
					tmp1.push_back(j);
				}
				else if (Array_metr[i][j] == Array_metr[j][i] && Array_metr[i][j] == 0) {

				}
				else {
					flag = false;
					break;
				}

			}
			if (flag == true) {
				if (vector_of_Q1.empty()) {
					//vector_of_Q1.push_back(names_of_row);
					vector_of_Q1.push_back(tmp1);
				}

				/*if (std::find(vector_of_Q1.begin(), vector_of_Q1.end(), names_of_row) == vector_of_Q1.end())
					vector_of_Q1.push_back(names_of_row);*/
				if (std::find(vector_of_Q1.begin(), vector_of_Q1.end(), tmp1) == vector_of_Q1.end())
					vector_of_Q1.push_back(tmp1);
				//cout << names_of_row << endl;
				//vector_of_Q1.push_back("1");
				/*names_of_row.clear();*/
				tmp1.clear();
			}

		}
		/*for (size_t n = 0; n < vector_of_Q1.size(); n++)
			cout << vector_of_Q1[n] << endl;*/

		cout << "___" << endl << endl;
		for (int i = 0; i < vector_of_Q1.size(); i++) {
			for (int j = 0; j < vector_of_Q1[i].size(); j++) {
				cout << vector_of_Q1[i][j] << "  ";
			}
			cout << endl;
		}

		vector<vector<int>> vector_of_Q2;
		vector<int> tmp2;
		for (int i = 0; i < m; i++) {
			bool flag = true;
			/*string names_of_row;*/
			for (int j = 0; j < n; j++) {
				if (Q2.Array_metr[i][j] == Q2.Array_metr[j][i] && Q2.Array_metr[i][j] == 1) {
					flag = true;
					/*	names_of_row += to_string(j);
							names_of_row += " ";*/
					tmp2.push_back(j);
				}
				else if (Q2.Array_metr[i][j] == Q2.Array_metr[j][i] && Q2.Array_metr[i][j] == 0) {

				}
				else {
					flag = false;
					break;
				}

			}
			if (flag == true) {
				if (vector_of_Q2.empty()) {
					//vector_of_Q1.push_back(names_of_row);
					vector_of_Q2.push_back(tmp2);
				}

				/*if (std::find(vector_of_Q1.begin(), vector_of_Q1.end(), names_of_row) == vector_of_Q1.end())
					vector_of_Q1.push_back(names_of_row);*/
				if (std::find(vector_of_Q2.begin(), vector_of_Q2.end(), tmp2) == vector_of_Q2.end())
					vector_of_Q2.push_back(tmp2);
				//cout << names_of_row << endl;
				//vector_of_Q1.push_back("1");
				/*names_of_row.clear();*/
				tmp2.clear();
			}

		}
		cout << "___" << endl << endl;
		for (int i = 0; i < vector_of_Q2.size(); i++) {
			for (int j = 0; j < vector_of_Q2[i].size(); j++) {
				cout << vector_of_Q2[i][j] << "  ";
			}
			cout << endl;
		}
		//TODO
		vector<vector<int>> vector_of_Q;



		for (int i = 0; i < vector_of_Q1.size(); i++) {

			int size_of_found = 0;
			for (int j = 0; j < vector_of_Q2.size(); j++) {
				bool smth_found = false;
				int size_of_found = 0;
				int size_in_begining = 0;
				vector<int> vector_of_clases;
				vector<int> tmp2 = vector_of_Q1[i];
				vector<int> tmp12 = vector_of_Q1[i];
				while (!tmp12.empty()) {
					size_in_begining = tmp2.size();
					while (!tmp2.empty()) {
						//cout << "here" << endl;
						//display(tmp2);
						//display(vector_of_Q2[j]);
						auto res = search(begin(vector_of_Q2[j]), end(vector_of_Q2[j]), begin(tmp2), end(tmp2));
						auto found = res != end(vector_of_Q2[j]);
						//cout << boolalpha << found;
						if (found == true && std::find(vector_of_Q.begin(), vector_of_Q.end(), tmp2) == vector_of_Q.end()) {//doesent work || i dont know why //+if not exist TODO
							vector_of_Q.push_back(tmp2);
							//cout << "+found" << endl;
							//display(tmp2);
							size_of_found += tmp2.size();
							smth_found = true;
							//cout << "size " << size_of_found << endl;

						}
						else { //
							tmp2.pop_back();//delete last and repeat
						}


					}

					tmp12.erase(tmp12.begin());//delete first & repeat
					tmp2 = tmp12;
				}
				//cout << "ENDED" << endl;


			}size_of_found = 0;

		}

		cout << "!!!!" << endl;

		for (int i = 0; i < vector_of_Q.size(); i++) {
			for (int j = 0; j < vector_of_Q[i].size(); j++) {
				cout << vector_of_Q[i][j] << "  ";

			}
			cout << endl;
		}

		vector<vector<int>> new_vector_of_Q;



		//for (int i = 0; i < vector_of_Q.size(); i++) {
		//	bool smd = false;
		//	for (int j = 0; j < vector_of_Q.size(); j++) {
		//		auto res = search(begin(vector_of_Q[i]), end(vector_of_Q[i]), begin(vector_of_Q[j]), end(vector_of_Q[j]));
		//		auto found = res != end(vector_of_Q2[j]);
		//		if (found == true) {//if exist
		//			smd = true;
		//			//vector_of_Q.erase(vector_of_Q[j].begin());
		//		}
		//		if (smd == false) {
		//			new_vector_of_Q.push_back(vector_of_Q[i]);
		//		}
		//	}
		//}


		cout << "?" << endl << endl;
		for (int h = 0; h < vector_of_Q.size(); h++) {
			for (int i = 0; i < vector_of_Q.size(); i++) {
				bool res_f = false;
				for (int j = 0; j < vector_of_Q.size(); j++) {
					//cout << "display 2" << endl;
					//display(vector_of_Q[i]);
					//display(vector_of_Q[j]);
					auto res2 = search(begin(vector_of_Q[i]), end(vector_of_Q[i]), begin(vector_of_Q[j]), end(vector_of_Q[j]));
					auto found22 = res2 != end(vector_of_Q[i]);
					if (i != j && found22 == true) {
						//cout << "oops" << endl;
						vector_of_Q.erase(vector_of_Q.begin() + j);
						res_f = true;

					}
				}
				if (res_f == false) {
					//cout << "push" << endl;
					//display(vector_of_Q[i]);
					new_vector_of_Q.push_back(vector_of_Q[i]);
				}
			}
		}



		cout << "______" << endl;
		for (int i = 0; i < vector_of_Q.size(); i++) {
			for (int j = 0; j < vector_of_Q[i].size(); j++) {
				cout << vector_of_Q[i][j] << "  ";

			}
			cout << endl;
		}

		vector<int> numbers_Q1;
		for (int i = 0; i < vector_of_Q1.size(); i++) {
			int num = 0;
			for (int j = 0; j < vector_of_Q.size(); j++) {

				bool is = false;
				for (int i_1 = 0; i_1 < vector_of_Q1[i].size(); i_1++) {

					for (int j_1 = 0; j_1 < vector_of_Q[j].size(); j_1++) {
						if (vector_of_Q1[i][i_1] == vector_of_Q[j][j_1]) {
							is = true;
						}



					}

					//is = false;
				}
				if (is == true) {
					num++;
				}

			}
			numbers_Q1.push_back(num);
			num = 0;
		}

		//	cout << "main" << endl;
			//display(numbers_Q1);


		vector<int> numbers_Q2;
		for (int i = 0; i < vector_of_Q2.size(); i++) {
			int num = 0;
			for (int j = 0; j < vector_of_Q.size(); j++) {

				bool is = false;
				for (int i_1 = 0; i_1 < vector_of_Q2[i].size(); i_1++) {

					for (int j_1 = 0; j_1 < vector_of_Q[j].size(); j_1++) {
						if (vector_of_Q2[i][i_1] == vector_of_Q[j][j_1]) {
							is = true;
						}



					}

					//is = false;
				}
				if (is == true) {
					num++;
				}

			}
			numbers_Q2.push_back(num);
			num = 0;
		}

		//cout << "main" << endl;
		//display(numbers_Q2);


	//	cout << "mostly" << endl;
		vector<int> numbers_Q1_MAX = laba5_help_makeMAX(vector_of_Q1, vector_of_Q);
		//	display(numbers_Q1_MAX);
		vector<int> numbers_Q2_MAX = laba5_help_makeMAX(vector_of_Q2, vector_of_Q);
		//display(numbers_Q2_MAX);






		float sum = 0;


		for (int i = 0; i < vector_of_Q1.size(); i++) {
			sum += 2.0*(vector_of_Q1[i].size() - numbers_Q1_MAX[i]) - numbers_Q1[i] + 1.0;
		}
		for (int i = 0; i < vector_of_Q2.size(); i++) {
			sum += 2.0*(vector_of_Q2[i].size() - numbers_Q2_MAX[i]) - numbers_Q2[i] + 1.0;
		}
		cout << "distance" << sum << endl;
	}

	vector < int> laba5_help_makeMAX(vector<vector<int>> vector_of_Q1, vector<vector<int>> vector_of_Q) {
		vector<int> numbers_Q1_max;
		for (int i = 0; i < vector_of_Q1.size(); i++) {
			int max = 0;
			for (int j = 0; j < vector_of_Q.size(); j++) {

				vector<int> tmp2 = vector_of_Q1[i];
				vector<int> tmp12 = vector_of_Q1[i];

				while (!tmp12.empty()) {

					while (!tmp2.empty()) {
						//cout << "THERE" << endl;
						//display(tmp2);
						//display(vector_of_Q[j]);
						auto res = search(begin(vector_of_Q[j]), end(vector_of_Q[j]), begin(tmp2), end(tmp2));
						auto found = res != end(vector_of_Q[j]);
						//cout << boolalpha << found << endl;
						if (found == true) {//doesent work || i dont know why //+if not exist TODO
							//cout << "+" << endl;
							//display(tmp2);
							if (tmp2.size() > max) {
								//cout << "max= " << tmp2.size() << endl;
								max = tmp2.size();
							}
							tmp2.pop_back();



						}
						else { //
							tmp2.pop_back();//delete last and repeat
						}


					}

					tmp12.erase(tmp12.begin());//delete first & repeat
					tmp2 = tmp12;
				}

				//cout << "ENDED" << endl;


			}
			numbers_Q1_max.push_back(max);

		}
		return numbers_Q1_max;
	}



	void find_measure_linear(relation_array_metr P2) {
		relation_array_metr P1_metr(this->m, this->n);
		relation_array_metr P2_metr(P2.m, P2.n);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				P1_metr.Array_metr[i][j] = this->Array_metr[i][j] - this->Array_metr[j][i];
				P2_metr.Array_metr[i][j] = P2.Array_metr[i][j] - P2.Array_metr[j][i];
			}
		}
		P1_metr.show();
		P2_metr.show();

		float distance = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				distance += 0.5*(fabs(P1_metr.Array_metr[i][j] - P2_metr.Array_metr[i][j]));
			}
		}
		cout << "distance= " << distance << endl;
	}
	bool help_check2() {

		bool match = true;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (Array_metr[i][j] == Array_metr[j][i] || (Array_metr[i][j] == 0 && Array_metr[j][i] == 0)) {
					//	cout << Array_metr[i][j] << " " << Array_metr[j][i] << " " << -1 * Array_metr[j][i] << endl;
					match = true;
				}
				else {
					match = false;
					return false;
					break;
				}

			}
		}
		return match;
	}

	bool is_Agreed() {
		bool agreed = false;
		bool match = help_check();
		bool match2 = help_check2();
		string type_str;
		if (match == true || match2 == true) {
			agreed = true;
		}

		return agreed;


	}
	void find_measure_metr(relation_array_metr P2) {
		double sum = 0;
		float num_W = 0;
		if (this->is_Agreed() && P2.is_Agreed()) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					if (this->Array_metr[i][j] != 100 && P2.Array_metr[i][j] != 100) {
						sum += fabs(this->Array_metr[i][j] - P2.Array_metr[i][j]);
					}
					else if (this->Array_metr[i][j] > P2.Array_metr[i][j]) {
						num_W++;
					}
					else {
						num_W++;
					}


				}
			}
			sum *= 0.5;
			cout << "distance= " << sum << " + " << num_W / 2 << "W" << endl;
		}
		else {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					if (this->Array_metr[i][j] != 100 && P2.Array_metr[i][j] != 100)
					{
						sum += fabs(this->Array_metr[i][j] - P2.Array_metr[i][j]);
					}
					else if (this->Array_metr[i][j] == 100 && P2.Array_metr[i][j] == 100) {
						sum += 0;
					}
					else if (this->Array_metr[i][j] > P2.Array_metr[i][j]) {
						num_W++;
					}
					else {
						num_W++;
					}
					/*else if (this->Array_metr[i][j] > P2.Array_metr[i][j]) {
						sum += this->Array_metr[i][j];
					}
					else {
						sum += P2.Array_metr[i][j];
					}*/
					if (this->Array_metr[j][i] != 100 && P2.Array_metr[j][i] != 100)
					{
						sum += fabs(this->Array_metr[j][i] - P2.Array_metr[j][i]);
					}
					else if (this->Array_metr[j][i] == 100 && P2.Array_metr[j][i] == 100) {
						sum += 0;
					}
					else if (this->Array_metr[j][i] > P2.Array_metr[j][i]) {
						num_W++;;
					}
					else {
						num_W++;
					}



				}
			}
			sum *= 0.25;
			cout << "distance= " << sum << " + " << num_W / 4 << "W" << endl;
		}

	}









	float **Array_metr;
};


//****************************************//
//lab4
//****************************************//



class mehanizm_of_choice {
public:
	virtual void find_solution() = 0;
	void show() {
		for (int i = 0; i < plural_alternative; i++) {
			for (int j = 0; j < plural_criteria; j++) {
				cout << plural_of_criterion_evaluations_of_alternatives[i][j] << " ";
			}cout << endl;
		}
		cout << endl;
	}

	//множина альтернатив,множина критеріїв, множина критерійних оцінок альтернатив

	int plural_alternative;
	int plural_criteria;
	float **plural_of_criterion_evaluations_of_alternatives;

};


class Paretto :public mehanizm_of_choice {
public:

	Paretto(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i<<"  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Paretto(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "Paretto" << endl;
		int m = plural_alternative, n = plural_criteria;
		relation_array ra(plural_alternative, plural_alternative);
		int count = -1;
		int count2 = 0;
		for (int i = 0; i < m; i++) {
			for (int k = 0; k < m; k++) {
				for (int j = 0; j < n; j++) {
					if (plural_of_criterion_evaluations_of_alternatives[i][j] >= plural_of_criterion_evaluations_of_alternatives[k][j]) {
						count = 0;
						if (plural_of_criterion_evaluations_of_alternatives[i][j] > plural_of_criterion_evaluations_of_alternatives[k][j]) {
							count2++;
						}
					}
					else {
						count = -1;
						break;
					}

				}
				if (count2 > 0 && count != -1) {
					ra.Array[i][k] = 1;
				}
				else {
					ra.Array[i][k] = 0;
				}
				count = -1;
				count2 = 0;

			}
		}

		ra.show();
		ra.find_majorant();
	}

};





class Sleyter :public mehanizm_of_choice {
public:

	Sleyter(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Sleyter(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "Sleyter" << endl;
		int m = plural_alternative, n = plural_criteria;
		relation_array ra(plural_alternative, plural_alternative);
		int count = -1;
		for (int i = 0; i < m; i++) {
			for (int k = 0; k < m; k++) {
				for (int j = 0; j < n; j++) {
					if (plural_of_criterion_evaluations_of_alternatives[i][j] > plural_of_criterion_evaluations_of_alternatives[k][j]) {
						count++;

					}
					else {
						count = -1;
						break;
					}

				}
				if (count >= 0) {
					ra.Array[i][k] = 1;
				}
				else {
					ra.Array[i][k] = 0;
				}
				count = -1;

			}
		}

		relation_array E(plural_alternative, plural_alternative, "diagonal");
		ra.relation_union(E);
		ra.show();
		ra.find_maximum();
	}

};






class The_best_result :public mehanizm_of_choice {
public:

	The_best_result(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	The_best_result(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "The Best result" << endl;
		int m = plural_alternative, n = plural_criteria;
		relation_array ra(plural_alternative, plural_alternative);
		//finding maximum
		float max = -1000;
		relation_array array_of_max(plural_alternative, 1);
		for (int i = 0; i < m; i++) {
			max = plural_of_criterion_evaluations_of_alternatives[i][0];
			for (int j = 0; j < n; j++) {
				if (plural_of_criterion_evaluations_of_alternatives[i][j] > max) {
					max = plural_of_criterion_evaluations_of_alternatives[i][j];
				}
			}
			array_of_max.Array[i][0] = max;
		}


		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				if (array_of_max.Array[i][0] >= array_of_max.Array[j][0]) {
					ra.Array[i][j] = 1;
				}
				else {
					ra.Array[i][j] = 0;
				}
			}
		}

		relation_array E(plural_alternative, plural_alternative, "diagonal");
		ra.relation_union(E);
		ra.show();
		ra.find_maximum();
	}

};




class Garanted_result :public mehanizm_of_choice {
public:

	Garanted_result(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Garanted_result(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "Garanted result" << endl;
		int m = plural_alternative, n = plural_criteria;
		relation_array ra(plural_alternative, plural_alternative);
		//finding minimum
		float min = 1000;
		relation_array array_of_min(plural_alternative, 1);
		for (int i = 0; i < m; i++) {
			min = plural_of_criterion_evaluations_of_alternatives[i][0];
			for (int j = 0; j < n; j++) {
				if (plural_of_criterion_evaluations_of_alternatives[i][j] < min) {
					min = plural_of_criterion_evaluations_of_alternatives[i][j];
				}
			}
			array_of_min.Array[i][0] = min;
		}
		//cout min
		// array_of_min.show();




		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				if (array_of_min.Array[i][0] >= array_of_min.Array[j][0]) {
					ra.Array[i][j] = 1;
				}
				else {
					ra.Array[i][j] = 0;
				}
			}
		}

		relation_array E(plural_alternative, plural_alternative, "diagonal");
		ra.relation_union(E);
		ra.show();
		ra.find_maximum();
	}

};



class Gurvic :public mehanizm_of_choice {
public:

	Gurvic(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Gurvic(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "Gurvic" << endl;
		float alfa;
		cout << "Input alfa:";
		cin >> alfa;
		int m = plural_alternative, n = plural_criteria;
		relation_array ra(plural_alternative, plural_alternative);
		//finding minimum
		float min = 1000;
		relation_array array_of_min(plural_alternative, 1);
		for (int i = 0; i < m; i++) {
			min = plural_of_criterion_evaluations_of_alternatives[i][0];
			for (int j = 0; j < n; j++) {
				if (plural_of_criterion_evaluations_of_alternatives[i][j] < min) {
					min = plural_of_criterion_evaluations_of_alternatives[i][j];
				}
			}
			array_of_min.Array[i][0] = min;
		}

		//array_of_min.show();

		//finding maximum
		float max = -1000;
		relation_array array_of_max(plural_alternative, 1);
		for (int i = 0; i < m; i++) {
			max = plural_of_criterion_evaluations_of_alternatives[i][0];
			for (int j = 0; j < n; j++) {
				if (plural_of_criterion_evaluations_of_alternatives[i][j] > max) {
					max = plural_of_criterion_evaluations_of_alternatives[i][j];
				}
			}
			array_of_max.Array[i][0] = max;
		}
		//cout max
		//array_of_max.show();



		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				if (alfa*array_of_min.Array[i][0] + (1 - alfa)*array_of_max.Array[i][0] > alfa * array_of_min.Array[j][0] + (1 - alfa)*array_of_max.Array[j][0]) {
					ra.Array[i][j] = 1;
				}
				else {
					ra.Array[i][j] = 0;
				}
			}
		}

		relation_array E(plural_alternative, plural_alternative, "diagonal");
		ra.relation_union(E);
		ra.show();
		ra.find_maximum();
	}

};




class Etalon :public mehanizm_of_choice {
public:

	Etalon(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Etalon(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "Etalon" << endl;
		float* etalon = new float[plural_criteria];
		cout << "Input Etalon: ";
		for (int i = 0; i < plural_criteria; i++) {
			cout << "i= ";
			cin >> etalon[i];
		}
		int m = plural_alternative, n = plural_criteria;
		relation_array ra(plural_alternative, plural_alternative);
		//finding minimum
		float min = 1000;
		relation_array array_of_etalon(plural_alternative, 1);
		for (int i = 0; i < m; i++) {
			array_of_etalon.Array[i][0] = 0;
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				//cout << fabs(plural_of_criterion_evaluations_of_alternatives[i][j] - etalon[j]) << endl;
				array_of_etalon.Array[i][0] += fabs(plural_of_criterion_evaluations_of_alternatives[i][j] - etalon[j]);
			}
			/*cout<< array_of_etalon.Array[i][0]
			array_of_etalon.Array[i][0] = 0;*/
		}
		//cout max
		array_of_etalon.show();



		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				if (array_of_etalon.Array[i][0] <= array_of_etalon.Array[j][0]) {
					ra.Array[i][j] = 1;
				}
				else {
					ra.Array[i][j] = 0;
				}
			}
		}

		/*relation_array E(plural_alternative, plural_alternative, "diagonal");
		ra.relation_union(E);*/
		ra.show();
		ra.find_maximum();
	}

};




//
//class Etalon :public mehanizm_of_choice {
//public:
//
//	Etalon(int amount_of_alter, int amount_of_criteria) {
//		float **matrix;
//		plural_alternative = amount_of_alter;
//		plural_criteria = amount_of_criteria;
//		int m = amount_of_alter, n = amount_of_criteria;
//		matrix = new  float*[amount_of_alter];
//		for (int z = 0; z < amount_of_alter; z++)
//			matrix[z] = new float[amount_of_criteria];
//		for (int i = 0; i < m; i++) {
//			//cout << "Input A" << i << "  ( ";
//			for (int j = 0; j < n; j++) {
//				cin >> matrix[i][j];
//			}
//			//cout << "  ) ";
//		}
//		plural_of_criterion_evaluations_of_alternatives = matrix;
//		show();
//	}
//
//
//
//	void find_solution() {
//		float* etalon = new float[plural_criteria];
//		cout << "Input Etalon: ";
//		for (int i = 0; i < plural_criteria; i++) {
//			cout << "i= ";
//			cin >> etalon[i];
//		}
//		int m = plural_alternative, n = plural_criteria;
//		relation_array ra(plural_alternative, plural_alternative);
//		//finding minimum
//		float min = 1000;
//		relation_array array_of_etalon(plural_alternative, 1);
//		for (int i = 0; i < m; i++) {
//			array_of_etalon.Array[i][0] = 0;
//		}
//		for (int i = 0; i < m; i++) {
//			for (int j = 0; j < n; j++) {
//				//cout << fabs(plural_of_criterion_evaluations_of_alternatives[i][j] - etalon[j]) << endl;
//				array_of_etalon.Array[i][0] += fabs(plural_of_criterion_evaluations_of_alternatives[i][j] - etalon[j]);
//			}
//			/*cout<< array_of_etalon.Array[i][0]
//			array_of_etalon.Array[i][0] = 0;*/
//		}
//		//cout max
//		array_of_etalon.show();
//
//
//
//		for (int i = 0; i < m; i++) {
//			for (int j = 0; j < m; j++) {
//				if (array_of_etalon.Array[i][0] <= array_of_etalon.Array[j][0]) {
//					ra.Array[i][j] = 1;
//				}
//				else {
//					ra.Array[i][j] = 0;
//				}
//			}
//		}
//
//		/*relation_array E(plural_alternative, plural_alternative, "diagonal");
//		ra.relation_union(E);*/
//		ra.show();
//	}
//
//};




class Convolution_criteria :public mehanizm_of_choice {
public:

	Convolution_criteria(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}

	Convolution_criteria(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}

	void find_solution() {
		cout << "Convolution criteria" << endl;
		float* vector_of_importace = new float[plural_criteria];

		float sum = 0;
		cout << "Input vector of importace of criteria : ";
		for (int i = 0; i < plural_criteria; i++) {
			cout << "i= ";
			cin >> vector_of_importace[i];
			sum += vector_of_importace[i];
		}

		if (sum != 1) {
			cout << "error --> sum must be 1" << endl;
		}
		else {
			int m = plural_alternative, n = plural_criteria;
			relation_array ra(plural_alternative, plural_alternative);
			relation_array_metr array_of_conv(plural_alternative, 1);
			float* array_of_conv2 = new float[plural_alternative];
			for (int i = 0; i < m; i++) {
				array_of_conv.Array[i][0] = 0.0;
				array_of_conv2[i] = 0.0;
			}
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					array_of_conv.Array[i][0] += plural_of_criterion_evaluations_of_alternatives[i][j] * vector_of_importace[j];
					array_of_conv2[i] += plural_of_criterion_evaluations_of_alternatives[i][j] * vector_of_importace[j];

				}

				/*cout << array_of_conv.Array[i][0] << endl;*/
				//array_of_etalon.Array[i][0] = 0;*/
			}
			//cout max



			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) {
					if (array_of_conv2[i] >= array_of_conv2[j]) {
						ra.Array[i][j] = 1;
					}
					else {
						ra.Array[i][j] = 0;
					}
				}
			}

			/*relation_array E(plural_alternative, plural_alternative, "diagonal");
			ra.relation_union(E);*/
			ra.show();
			ra.find_maximum();
		}
	}
};





class Lexicographical :public mehanizm_of_choice {
public:

	Lexicographical(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Lexicographical(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		cout << "lexicographical" << endl;
		int* order_of_importace = new int[plural_criteria];
		cout << "Input order of importace of criteria : ";
		for (int i = 0; i < plural_criteria; i++) {
			cout << "Q= ";
			cin >> order_of_importace[i];

		}
		int m = plural_alternative, n = plural_criteria;

		relation_array ra(plural_alternative, plural_alternative, "empty");
		for (int i = 0; i < m; i++) {
			for (int i1 = 0; i1 < m; i1++) {
				for (int j = 0; j < n; j++) {
					if (plural_of_criterion_evaluations_of_alternatives[i][order_of_importace[j]] > plural_of_criterion_evaluations_of_alternatives[i1][order_of_importace[j]]) {
						ra.Array[i][i1] = 1;
						break;
					}
					if (plural_of_criterion_evaluations_of_alternatives[i][order_of_importace[j]] < plural_of_criterion_evaluations_of_alternatives[i1][order_of_importace[j]]) {
						ra.Array[i][i1] = 0;
						break;
					}

				}
			}
		}


		/*relation_array E(plural_alternative, plural_alternative, "diagonal");
		ra.relation_union(E);*/
		ra.show();
		ra.find_majorant();





	}

};






class Main_criteria :public mehanizm_of_choice {
public:

	Main_criteria(int amount_of_alter, int amount_of_criteria) {
		float **matrix;
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		int m = amount_of_alter, n = amount_of_criteria;
		matrix = new  float*[amount_of_alter];
		for (int z = 0; z < amount_of_alter; z++)
			matrix[z] = new float[amount_of_criteria];
		for (int i = 0; i < m; i++) {
			//cout << "Input A" << i << "  ( ";
			for (int j = 0; j < n; j++) {
				cin >> matrix[i][j];
			}
			//cout << "  ) ";
		}
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}
	Main_criteria(int amount_of_alter, int amount_of_criteria, float **matrix) {
		plural_alternative = amount_of_alter;
		plural_criteria = amount_of_criteria;
		plural_of_criterion_evaluations_of_alternatives = matrix;
		show();
	}


	void find_solution() {
		int* order_of_importace = new int[plural_criteria];
		int position = -1;
		cout << "Method main criteria \n Input criteria|| main criteria =99 : ";
		for (int i = 0; i < plural_criteria; i++) {
			cout << "Q" << i << "=";
			cin >> order_of_importace[i];
			if (order_of_importace[i] == 99) {
				position = i;
			}
		}


		int m = plural_alternative, n = plural_criteria;

		relation_array ra(plural_alternative, plural_alternative, "empty");
		for (int i = 0; i < m; i++) {
			for (int i1 = 0; i1 < m; i1++) {
				for (int j = 0; j < n; j++) {
					if (plural_of_criterion_evaluations_of_alternatives[i][position] >= plural_of_criterion_evaluations_of_alternatives[i1][position]) {
						bool is_true = false;
						for (int k = 0; k < plural_criteria; k++) {
							if (k != position) {
								if (plural_of_criterion_evaluations_of_alternatives[i][k] >= order_of_importace[k] && plural_of_criterion_evaluations_of_alternatives[i1][k] >= order_of_importace[k]) {
									is_true = true;
								}
								else {
									is_true = false;
									break;
								}


							}
							//cout << plural_of_criterion_evaluations_of_alternatives[i][k] << "  " << plural_of_criterion_evaluations_of_alternatives[i1][k] << "  " << order_of_importace[k] << endl;

						}
						if (is_true == true) {

							ra.Array[i][i1] = 1;
						}
						else {
							ra.Array[i][i1] = 0;
						}
					}
					if (plural_of_criterion_evaluations_of_alternatives[i][position] < plural_of_criterion_evaluations_of_alternatives[i1][position]) {
						ra.Array[i][i1] = 0;
						//break;
					}

				}
			}
		}

		ra.show();
		vector<int> the_same;
		for (int i = 0; i < m; i++) {
			bool row = true;
			for (int j = 0; j < m; j++) {
				if (ra.Array[i][j] == 0) {
					row = false;

				}
				else {
					row = true;
					break;
				}
			}
			if (row == false) {
				the_same.push_back(i);
			}
		}



		ra.relation_narrowing(the_same);
		ra.show();
		ra.find_maximum();





	}

};


// TODO



class group_relation : public virtual relation_array_metr {
public:
	//множина альтернатив, множина принципів узгодження, множина експертних оцінок альтернатив.
	float *plural_alternative;
	relation_array_metr expert_marks;
	virtual void find_solution() = 0;
	void show() {
		for (int i = 0; i < expert_marks.m; i++) {
			for (int j = 0; j < expert_marks.n; j++) {
				cout << expert_marks.Array_metr[i][j] << " ";
			}cout << endl;
		}
		cout << endl;
	}
};

class majority_votes : public group_relation {
public:
	vector<relation_array_metr> P;
	majority_votes(vector<relation_array_metr> P1) {
		P = P1;
	}

	void find_solution() {

		relation_array_metr R(P[0].m, P[0].n, "empty");
		for (int k = 0; k < P.size(); k++) {
			int *sum = new int[P[0].m];
			for (int i = 0; i < P[0].n; i++) {
				sum[i] = 0;
			}
			for (int i = 0; i < P[0].m; i++) {
				for (int j = 0; j < P[0].n; j++) {
					if (P[k].Array_metr[i][j] == 1) {
						sum[i]++;
					}
					//R[k].Array_metr[i][j] = P[k].Array_metr[i][j] - P[k].Array_metr[j][i];
				}
			}

			for (int i = 0; i < P[0].n; i++) {
				if (sum[i] == 4) {
					R.Array_metr[0][k] = i + 1;
				}
				if (sum[i] == 3) {
					R.Array_metr[1][k] = i + 1;
				}
				if (sum[i] == 2) {
					R.Array_metr[2][k] = i + 1;
				}
				if (sum[i] == 1) {
					R.Array_metr[3][k] = i + 1;
				}
				if (sum[i] == 0) {
					R.Array_metr[4][k] = i + 1;
				}

			}


		}
		cout << "!!!" << endl;
		R.show();



		int *kilkist_povtorov = new int[R.m];
		for (int j = 0; j < R.m; j++) {
			kilkist_povtorov[j] = count(R.Array_metr[0], R.Array_metr[0] + R.m, R.Array_metr[0][j]);
		}
		for (int j = 0; j < R.m; j++) {
			cout << "num: " << kilkist_povtorov[j] << endl;
		}
		int max = kilkist_povtorov[0];
		int position = 0;
		for (int j = 0; j < R.m; j++) {
			if (kilkist_povtorov[j] > max) {
				max = kilkist_povtorov[j];
				position = j;
			}
		}
		cout << "Kondorse max alternative is x" << R.Array_metr[0][position] << endl;

	}
};






class Borda : public group_relation {
public:
	vector<relation_array_metr> P;
	/*Borda(int m, int n, float **arr) {
		relation_array_metr expert_marks2(m, n, arr);
		expert_marks = expert_marks2;
	}*/
	Borda(vector<relation_array_metr> P1) {
		P = P1;
	}

	void find_solution() {


		relation_array_metr R(P[0].m, P[0].n, "empty");
		for (int k = 0; k < P.size(); k++) {
			int *sum = new int[P[0].m];
			for (int i = 0; i < P[0].n; i++) {
				sum[i] = 0;
			}
			for (int i = 0; i < P[0].m; i++) {
				for (int j = 0; j < P[0].n; j++) {
					if (P[k].Array_metr[i][j] == 1) {
						sum[i]++;
					}
					//R[k].Array_metr[i][j] = P[k].Array_metr[i][j] - P[k].Array_metr[j][i];
				}
			}

			for (int i = 0; i < P[0].n; i++) {
				if (sum[i] == 4) {
					R.Array_metr[0][k] = i + 1;
				}
				if (sum[i] == 3) {
					R.Array_metr[1][k] = i + 1;
				}
				if (sum[i] == 2) {
					R.Array_metr[2][k] = i + 1;
				}
				if (sum[i] == 1) {
					R.Array_metr[3][k] = i + 1;
				}
				if (sum[i] == 0) {
					R.Array_metr[4][k] = i + 1;
				}

			}


		}
		cout << "!!!" << endl;
		R.show();


		int *sum_borda = new int[R.m];
		for (int i = 0; i < R.n; i++) {
			sum_borda[i] = 0;
		}
		for (int i = 0; i < R.m; i++) {
			for (int j = 0; j < R.n; j++) {
				sum_borda[int(R.Array_metr[i][j]) - 1] += R.m - (i + 1);

			}
		}
		for (int i = 0; i < R.n; i++) {
			cout << sum_borda[i] << endl;
		}
		int max = sum_borda[0];
		int position = 0;
		for (int j = 0; j < R.m; j++) {
			if (sum_borda[j] > max) {
				max = sum_borda[j];
				position = j;
			}
		}
		cout << "Borda max alternative is x" << position + 1 << endl;
	}
};





class mediana_kemeni : public group_relation {
public:
	vector<relation_array_metr> P;
	mediana_kemeni(vector<relation_array_metr> P1) {
		P = P1;
	}

	void find_solution() {
		//change to metr
		int amount_of_alternative = 5;
		relation_array_metr p1(amount_of_alternative, amount_of_alternative, "empty");
		relation_array_metr p2(amount_of_alternative, amount_of_alternative, "empty");
		relation_array_metr p3(amount_of_alternative, amount_of_alternative, "empty");
		relation_array_metr p4(amount_of_alternative, amount_of_alternative, "empty");
		relation_array_metr p5(amount_of_alternative, amount_of_alternative, "empty");
		/*relation_array_metr p5(5, 5, "input");
		p5.show();*/
		vector<relation_array_metr> P_metr{ p1,p2,p3,p4,p5 };

		for (int k = 0; k < P.size(); k++) {
			for (int i = 0; i < P[0].m; i++) {
				for (int j = 0; j < P[0].n; j++) {
					P_metr[k].Array_metr[i][j] = P[k].Array_metr[i][j] - P[k].Array_metr[j][i];
				}
			}
		}
		cout << "metr" << endl;
		for (int k = 0; k < P_metr.size(); k++) {
			P_metr[k].show();
		}
		//find R

		relation_array_metr R(P[0].m, P[0].n, "empty");

		for (int i = 0; i < P[0].m; i++) {
			for (int j = 0; j < P[0].n; j++) {
				if (i == j) {
					R.Array_metr[i][j] = 0;
				}
				else {
					for (int k = 0; k < P_metr.size(); k++) {
						R.Array_metr[i][j] += 1 - P_metr[k].Array_metr[i][j];
					}
				}
			}
		}
		R.show();
		/*
		}*/
		//delete from R
		int num = 0;
		vector<int> digits_num;
		while (R.m != 0) {
			cout << "iteration " << num + 1 << endl;
			float *array_sum = new float[R.m];
			for (int i = 0; i < R.m; i++) {
				int sum = 0;
				for (int j = 0; j < R.n; j++) {
					sum += R.Array_metr[i][j];
				}
				array_sum[i] = sum;
			}
			for (int i = 0; i < R.m; i++) {
				cout << array_sum[i] << endl;
			}

			float min = array_sum[0];
			int pos = 0;
			for (int i = 0; i < R.m; i++) {
				if (array_sum[i] < min) {
					min = array_sum[i];
					pos = i;
				}
			}
			cout << "zuzh" << pos + 1 << endl;
			R.relation_narrowing({ pos });
			digits_num.push_back(pos + 1);
			R.show();
			num++;
			cout << "R.m " << R.m << " R.n " << R.n << endl;

		}
		for (int i = 0; i < digits_num.size(); i++) {
			cout << digits_num[i] << "  ";
		}
		for (int i = 0; i < digits_num.size(); i++) {
			for (int j = 0; j < i; j++) {
				if (digits_num[i] >= digits_num[j]) {
					digits_num[i]++;
				}
			}
		}
		cout << endl << "KEMENI ";
		for (int i = 0; i < digits_num.size(); i++) {
			cout << digits_num[i] << "  ";
		}
	}
};



//****************************************//
//****************************************//
void Bayes_Laplace(relation_array_metr mtx, relation_array_metr probab) {
	relation_array_metr a(mtx.m, 1, "empty");

	for (int i = 0; i < mtx.m; i++) {
		for (int j = 0; j < mtx.n; j++) {
			a.Array_metr[i][0] += mtx.Array_metr[i][j] * probab.Array_metr[0][j];
		}
	}

	a.show();
	int position = 0;
	float max = a.Array_metr[0][0];
	for (int i = 0; i < mtx.m; i++) {
		if (a.Array_metr[i][0] > max) {
			max = a.Array_metr[i][0];
			position = i;
		}
	}

	cout << "Bayes_Laplace max is alternative " << position + 1 << " with " << max << endl;
}

void Valda(relation_array_metr mtx) {

	relation_array_metr mins(mtx.m, 1, "empty");
	for (int i = 0; i < mtx.m; i++) {
		mins.Array_metr[i][0] = mtx.Array_metr[i][0];
		for (int j = 0; j < mtx.n; j++) {
			if (mtx.Array_metr[i][j] < mins.Array_metr[i][0]) {
				mins.Array_metr[i][0] = mtx.Array_metr[i][j];
			}
		}
	}
	mins.show();

	int position = 0;
	float max = mins.Array_metr[0][0];
	for (int i = 0; i < mtx.m; i++) {
		if (mins.Array_metr[i][0] > max) {
			max = mins.Array_metr[i][0];
			position = i;
		}
	}

	cout << "Valda max is alternative " << position + 1 << " with " << max << endl;
}
void Seavedge(relation_array_metr mtx) {
	mtx.show();
	relation_array_metr maxs(mtx.n, 1, "empty");
	for (int j = 0; j < mtx.n; j++) {

		maxs.Array_metr[j][0] = mtx.Array_metr[0][j];
		for (int i = 0; i < mtx.m; i++) {

			if (mtx.Array_metr[i][j] > maxs.Array_metr[j][0]) {
				maxs.Array_metr[j][0] = mtx.Array_metr[i][j];
			}
		}
	}
	maxs.show();
	relation_array_metr R(mtx.m, mtx.n, "empty");
	for (int i = 0; i < mtx.m; i++) {
		for (int j = 0; j < mtx.n; j++) {
			R.Array_metr[i][j] = maxs.Array_metr[j][0] - mtx.Array_metr[i][j];
		}
	}
	cout << "Risk matrix: " << endl;
	R.show();


	relation_array_metr maxs_r(R.m, 1, "empty");
	for (int i = 0; i < R.m; i++) {
		maxs_r.Array_metr[i][0] = R.Array_metr[i][0];
		for (int j = 0; j < R.n; j++) {
			if (R.Array_metr[i][j] > maxs_r.Array_metr[i][0]) {
				maxs_r.Array_metr[i][0] = R.Array_metr[i][j];
			}
		}
	}
	maxs_r.show();

	int position = 0;
	float min = maxs_r.Array_metr[0][0];
	for (int i = 0; i < maxs_r.m; i++) {
		if (maxs_r.Array_metr[i][0] < min) {
			min = maxs_r.Array_metr[i][0];
			position = i;
		}
	}

	cout << "Seavedge min(max(R)) is alternative " << position + 1 << " with " << min << endl;


}
void Gurvic(relation_array_metr mtx, float alpha) {
	relation_array_metr mins(mtx.m, 1, "empty");
	for (int i = 0; i < mtx.m; i++) {
		mins.Array_metr[i][0] = mtx.Array_metr[i][0];
		for (int j = 0; j < mtx.n; j++) {
			if (mtx.Array_metr[i][j] < mins.Array_metr[i][0]) {
				mins.Array_metr[i][0] = mtx.Array_metr[i][j];
			}
		}
	}
	mins.show();

	relation_array_metr maxs(mtx.m, 1, "empty");
	for (int i = 0; i < mtx.m; i++) {

		maxs.Array_metr[i][0] = mtx.Array_metr[i][0];
		for (int j = 0; j < mtx.n; j++) {

			if (mtx.Array_metr[i][j] > maxs.Array_metr[i][0]) {
				maxs.Array_metr[i][0] = mtx.Array_metr[i][j];
			}
		}
	}
	maxs.show();
	relation_array_metr s(mtx.m, 1, "empty");
	for (int i = 0; i < mtx.m; i++) {
		s.Array_metr[i][0] = alpha * mins.Array_metr[i][0] + (1 - alpha)*maxs.Array_metr[i][0];
	}
	s.show();

	int position = 0;
	float max = s.Array_metr[0][0];
	for (int i = 0; i < s.m; i++) {
		if (s.Array_metr[i][0] > max) {
			max = s.Array_metr[i][0];
			position = i;
		}
	}

	cout << "Gurvic max is alternative " << position + 1 << " with " << max << endl;


}




void lab7() {
	int V, P;
	cout << "input amount of strategies V: ";
	cin >> V;
	cout << "input amount of demand P(poput): ";
	cin >> P;
	cout << "input Matrix" << endl;
	relation_array_metr MATRIX(V, P, "input");
	MATRIX.show();
	cout << "input p(probability) " << endl;
	relation_array_metr p(1, P, "input");
	cout.setf(ios::fixed);
	Bayes_Laplace(MATRIX, p);
	Valda(MATRIX);
	Seavedge(MATRIX);
	float num = 0;
	while (num <= 1) {
		cout << "!!!num " << num << endl;
		Gurvic(MATRIX, num);
		num += 0.2;
	}
}


class A {
public:
	void foo() const {
		std::cout << "A";
	}
};
class B {
public:
	void foo() const {
		std::cout << "B";
	}
};
class C : public A, public B {
	using A::foo;
};
struct Foo {
	int x;
	int y;

};

int fn(int v) {
	if (v == 1 || v == 0) {
		return 1;
	}
	if (v % 2 == 0) {
		return fn(v / 2) + 2;
	}
	else return fn(v - 1) + 3;
}
int x = 20;

namespace outer {
	int x = 10;
	namespace inner {
		int z = x;
	}
}





int main() {

	cout << x + outer::x << endl;
	system("pause");
	return 0;
	/*MEMORYSTATUSEX memInfo;
	memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo);
	DWORDLONG totalVirtualMem = memInfo.ullTotalPageFile;
	DWORDLONG physMemUsed = memInfo.ullTotalPhys - memInfo.ullAvailPhys;
	cout <<"memory in the begining "<< physMemUsed << endl;

	relation_cat P(5, 5, "input");
	relation_cat Q(5, 5, "input");
	relation_cat R(5, 5, "input");

	P.show();
	Q.show();
	R.show();
	P.calculate(Q, R);
	cout << "resul :" << endl;
	P.show();

	MEMORYSTATUSEX memInfo2;
	memInfo2.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo2);
	DWORDLONG totalVirtualMem2 = memInfo2.ullTotalPageFile;

	DWORDLONG physMemUsed2 = memInfo2.ullTotalPhys - memInfo2.ullAvailPhys;

	cout << "memory in the end " << physMemUsed2 << endl;

	cout << "memory used " << (physMemUsed2 - physMemUsed )/1024<<"kb"<< endl;
	*/


	//relation_array_metr ram(5, 5, "input");
	//ram.show();
	//cout << ram.determining_type_of_relation() << endl;
	//cout << ram.check_matched() << endl;
	//ram.check_transitive();


	//relation_array_metr ram2(5, 5, "input");
	//ram2.show();
	//cout << ram2.determining_type_of_relation() << endl;
	//cout << ram2.check_matched() << endl;
	//ram2.check_transitive();


	//relation_array_metr ram3(5, 5);
	//ram3 = ram; //union 
	//relation_array_metr ram4(5, 5);
	//ram4 = ram; //intersection
	//relation_array_metr ram5(5, 5);
	//ram5 = ram; //difference
	//relation_array_metr ram6(5, 5);
	//ram6 = ram; //composition

	//ram3.relation_union(ram2);
	//cout << "relation union:" << endl;
	//ram3.show();

	//ram4.relation_intersection(ram2);
	//cout << "relation intersection:" << endl;
	//ram4.show();

	//ram5.relation_difference(ram2);
	//cout << "relation difference:" << endl;
	//ram5.show();

	//ram6.relation_composition(ram2);
	//cout << "relation composition:" << endl;
	//ram6.show();

	//



	//relation_array arr(5, 5, "input");
	//arr.show();

	//arr2 = arr.find_mutual_reach();
	//arr2.show();
	//relation_array arr3 = arr.factorization_by_the_given_relation_equivalence(arr2);

	//cout << "out" << endl << endl << endl;
	//arr3.show();

//	arr3.check_type();


	//relation_array_metr ram2(4, 4, "input");
	//ram2.show();
	//ram2.find_maximum();
	//ram2.find_minimum();
	//ram2.find_majorant();
	//ram2.find_minorant();



	/*float **matrix;
	int plural_alternative,	plural_criteria;

	cout << "Input plural alternative: " ;
	cin >> plural_alternative;
	cout << "Input plural criteria:  " ;
	cin >> plural_criteria;
	matrix = new  float*[plural_alternative];
	for (int z = 0; z < plural_alternative; z++)
		matrix[z] = new float[plural_criteria];
	cout << "Input matrix" << endl;
	for (int i = 0; i < plural_alternative; i++) {
		for (int j = 0; j < plural_criteria; j++) {
			cin >> matrix[i][j];
		}
	}*/
	////4laba
	//Paretto pareto(plural_alternative, plural_criteria, matrix);
	//pareto.find_solution();
	//Sleyter sleyter(plural_alternative, plural_criteria, matrix);
	//sleyter.find_solution();
	//The_best_result the_best_result(plural_alternative, plural_criteria, matrix);
	//the_best_result.find_solution();
	//Garanted_result garanted_res(plural_alternative, plural_criteria, matrix);
	//garanted_res.find_solution();
	//Gurvic gurvic(plural_alternative, plural_criteria, matrix);
	//gurvic.find_solution();
	//Etalon etalon(plural_alternative, plural_criteria, matrix);
	//etalon.find_solution();
	//Convolution_criteria convolution_criteria(plural_alternative, plural_criteria, matrix);
	//convolution_criteria.find_solution();
	//Lexicographical lexicographical(plural_alternative, plural_criteria, matrix);
	//lexicographical.find_solution();
	//Main_criteria main_criteria(plural_alternative, plural_criteria, matrix);
	//main_criteria.find_solution();
	////end of 4 laba


	//vector<string> v1 = { "A","B","C" };
	//vector<string> v2 = { "X","Y","A","B","C","D" };

	//auto res = search(begin(v2), end(v2), begin(v1), end(v1));
	//cout << *res << endl;

	//if (std::find(proba.begin(), proba.end(), proba1) != proba.end()) {
	//	cout << "found" << endl;
	//}
	//
	//for (std::vector<int>::const_iterator iter_rule_list = ruleList[0].begin(); iter_rule_list != ruleList[0].end(); ++iter_rule_list)
	//{
	//	int const value = *iter_rule_list;
	//	if (std::find(ntList.begin(), ntList.end(), value) == ntList.end())
	//	{
	//		// value missing in ntList
	//	}
	//}
	//

//equivavalence
	/*relation_array_metr Q1(5, 5, "input");
	Q1.show();
	relation_array_metr Q2(5, 5, "input");
	Q2.show();
	Q1.find_measure_of_proximity_equivalence(Q2);
*/



//int NN = 5;

//relation_array_metr Q1(NN, NN, "input");
//Q1.show();
//relation_array_metr Q2(NN, NN, "input");
//Q2.show();
////Q1.find_measure_linear(Q2);
//Q1.find_measure_metr(Q2);

//lab6
	//float **matrix;
	//int plural_alternatives=5,	plural_criteria=5;
	///*cout << "Input plural alternative: " ;
	//cin >> plural_alternatives;
	//cout << "Input plural criteria:  " ;
	//cin >> plural_criteria;*/
	//matrix = new  float*[plural_alternatives];
	//for (int z = 0; z < plural_alternatives; z++)
	//	matrix[z] = new float[plural_criteria];
	//cout << "Input matrix" << endl;
	//for (int i = 0; i < plural_alternatives; i++) {
	//	for (int j = 0; j < plural_criteria; j++) {
	//		cin >> matrix[i][j];
	//	}
	//}
	//majority_votes majority(plural_alternatives, plural_criteria, matrix);
	//majority.show();
	//majority.find_solution();

	//Borda borda(plural_alternatives, plural_criteria, matrix);
	//borda.show();
	//borda.find_solution();
//lab 6
//cout << "Input matrix" << endl;
//	int amount_of_alternative = 5;
//	relation_array_metr p1(amount_of_alternative, amount_of_alternative, "input");
//	relation_array_metr p2(amount_of_alternative, amount_of_alternative, "input");
//	relation_array_metr p3(amount_of_alternative, amount_of_alternative, "input");
//	relation_array_metr p4(amount_of_alternative, amount_of_alternative, "input");
//	relation_array_metr p5(amount_of_alternative, amount_of_alternative, "input");
//	p1.show();
//	p2.show();
//	p3.show();
//	p4.show();
//	p5.show();
//	
//	//p5.show();*/
//	vector<relation_array_metr> expert{ p1,p2,p3,p4,p5 };
//	//Borda borda(amount_of_alternative, amount_of_alternative,p1.Array_metr);
//	Borda borda(expert);
//	borda.find_solution();
//
//	majority_votes majority(expert);
//	majority.find_solution();
//
//	mediana_kemeni kemeni(expert);
//	kemeni.find_solution();
	malloc(sizeof(0));

	lab7();
	system("pause");
	return 0;

}