#pragma once
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

#include <vector>

// Вспомогательные функции:
double compute_tetra_volume();

void write_tecplot(mesh, data, fname, var_names, time = 0.0);

class Mesh{
//	пока можно всё делать в public
public:
	void read_starcd(const string& path, const int scale = 1);
	// Нужно придумать удобные типы для всех массивов, например
	vector<vector<int>> bcface_vert_lists;
	// ...

};
