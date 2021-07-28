#ifndef MESH_H_
#define MESH_H_

#include "header.h"

using namespace std;
/*
 * Header for Mesh class
 * Чтобы не напрягаться из-за ручной работы с памятью
 * нужно использовать для основных данных стандартные контейнеры:
 * vector, set, map и т.п. Важно помнить:
 * Элементы в set, map автоматически сортируются
 * можно создавать vector<vector<int>> а , тогда a[i] могут быть
 * разной длины - это удобно
 * Можно за один проход заполнять коллекции -
 * не нужно вычислять заранее, какие должны быть размеры
 *
 * Возможно в будущем нужно будет переделать класс так,
 * чтобы все структуры были в private,
 * а доступ к нужным элементам выполнялся через методы Get...
 *
 * Нужно иметь ввиду, что в будущем нужно будет всё обобщить на сетку
 * с разными типами ячеек (тэтраэдры, призмы и т.п.)
 *
 * construct faces of cell - желательно вынести в отдельный метод
 *
 * Декомпозицию на тэтраэдры -
 * желательно вынести в отдельную функцию (# 1st tetra и т.д.)
 */

class Mesh {
public:
	// Constructor
	Mesh(const string& path, double scale = 1.0);
	void read_starcd(const string& path, double scale = 1.0);

	double compute_tetra_volume(vector < vector < double > > tetra);
	vector < vector < int > > get_faces_for_cell(int ic);

	int nbf;
	int nbc;
	int nc;
	int nv;
	int nf;

	vector < vector < int > > bcface_vert_lists;
	vector < set < int > > bcface_vert_set_lists; // TODO: Don't need set
	vector < int > bcface_bctype;
	vector < vector < int > > vert_list_for_cell;
	vector < vector < int > > bf_for_each_bc;
	vector < vector < double > > vert_coo;
	vector < vector < double > > cell_center_coo;
	vector < double > cell_volumes;
	vector < vector < int > > cell_list_for_vertex;
	vector < int > cell_num_for_vertex;
	vector < vector < int > > cell_neighbors_list;
	vector < vector < int > > face_vert_list;
	vector < vector < int > > cell_face_list;
	vector < double > face_areas;
	vector < vector < double > > face_normals;
	vector < vector < int > > cell_face_normal_direction;
	vector < vector < int > > bound_face_info;
	vector < double > cell_diam;
	vector < vector < double > > face_centers;

	void write_tecplot(vector < vector <double> > data, string filename,
			vector <string> var_names, double time = 0.0);
};

#endif /* MESH_H_ */
