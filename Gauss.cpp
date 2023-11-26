#include <iostream>//заголовочный файл
#include <iomanip>
#include <stdlib.h>
#include <math.h>

using namespace std;

float determinant(float **matrix, int n);
float * sideDeters(float **matrix, float * b, int n);
float * methodGauss(float **matrix, int n, int m);
void display(float **matrix, int n, int m);
void mtrxCopy(float ** source, float ** target, int rows, int columns);

void FreeMem(float **matrix, int n);
void PrintMtx(float **matrix, int n);
void TransponMtx(float **matrix, float **tMatrix, int n);
void Get_matr(float **matrix, int n, float **temp_matrix, int indRow, int indCol);
int Det(float **matrix, int n);

int main()//начало выполения программы
{
	setlocale(LC_CTYPE, "rus");//перевод на русский
	int i, j, n, m;//объявление переменных
	double *x;
	cout << "Введите число уравнений: \n";//отображение
	cin >> n;//ввод данных отпользователя
	cout << "Введите число неизвестных: \n";//отображение
	cin >> m;//ввод данных от пользователя
	m += 1;
	x = new double[n];
	float **matrix = new float*[n];
	for (i = 0; i < n; i++)
		matrix[i] = new float[m];

	float **obr_matrix = new float *[n];
	for (int i = 0; i < n; i++)
		obr_matrix[i] = new float[n];

	float **tobr_matrix = new float *[n];
	for (int i = 0; i < n; i++) 
		tobr_matrix[i] = new float[n];

	//инициализация
	for (i = 0; i < n; i++) {
		for (j = 0; j < m-1; j++)
		{
			cout << "a" << "[" << i+1 << "]"<<"[" << j+1 << "]= ";//отображение
			cin >> matrix[i][j];
		}
		cout << "b" << "[" << i+1 << "]= ";//отображение
		cin >> matrix[i][m-1];
	}
	//выводим массив
	cout << "Исходная матрица: \t " << endl;
	display(matrix, n, m);//вызываем процедура
	cout << endl;

	//Вызов метода Гаусса
	float  *back = new float[m];
	back = methodGauss(matrix, n, m);//вызываем процедуру
	//Выводим решения

	for (i = 0; i < n; i++)
		//Выводим 4 знака после запятой (убрать, если нужна большая точность)
		cout << "Ответ для метода Гаусса" << "[" << i + 1 << "]= " << setiosflags(ios::fixed) << setprecision(4) <<back[i] << endl;
	cout << endl;
	int det;

	PrintMtx(matrix, n);
	det = Det(matrix, n);
	cout << "Определитель матрицы = " << det << endl;
	if (det) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int m = n - 1;
				float **temp_matrix = new float *[m];
				for (int k = 0; k < m; k++)
					temp_matrix[k] = new float[m];
				Get_matr(matrix, n, temp_matrix, i, j);
				obr_matrix[i][j] = pow(-1.0, i + j + 2) * Det(temp_matrix, m) / det;
				FreeMem(temp_matrix, m);
			}
		}
	}
	else
		cout << "Т.к. определитель матрицы = 0,\nто матрица вырожденная и обратной не имеет!!!)" << endl;
	//Транспонирование матрицы
	TransponMtx(obr_matrix, tobr_matrix, n);
	cout << "Обратная матрица: \t" << endl;
	
	cout<< "Основной детерминант = " << determinant(matrix, n) << endl;

	float * b = new float[n];

	for (i = 0; i < n; i++)
	{
		cout << "b" << "[" << i + 1 << "]= ";
		cin >> b[i];
	}
	
	float * deters = sideDeters(matrix, b, n);
	for (i = 0; i < n; i++) {
		cout << "\n Детерминант побочный" << i+ 1 << " = " << deters[i] << endl;
	}
	cout << endl;

	//Печать обратной матрицы после транспонирования
	PrintMtx(tobr_matrix, n);
	FreeMem(tobr_matrix, n);
	FreeMem(matrix, n);
	FreeMem(obr_matrix, n);

	delete[] matrix;//освобождение памяти,выделенной под массив
	cout << "Для продолжения нажмите любую клавишу ...";
	cin.get(); cin.get();//ожидание ввода пользователем любого символа,после считывания которого программа завершается
	return 0;
}
float determinant(float **matrix, int n)
{
	int i, j, k, r;
	double l, M, max, det = 1;

	float **matrix1 = new float*[n];
	mtrxCopy(matrix, matrix1, n, n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			matrix1[i][j] = matrix[i][j];
		}
	}

	for (k = 0; k<n; k++)
	{
		max = fabs(matrix1[k][k]);
		r = k;
		for (i = k + 1; i<n; i++)
			if (fabs(matrix1[i][k])>max)
			{
				max = fabs(matrix1[i][k]);
				r = i;
			}
		if (r != k)
		{
			det = -det;
		}
		for (j = 0; j<n; j++)
		{
			l = matrix1[k][j];
			matrix1[k][j] = matrix1[r][j];
			matrix1[r][j] = l;
		}
		for (i = k + 1; i<n; i++)
			for (M = matrix1[i][k] / matrix1[k][k], j = k; j<n; j++)
				matrix1[i][j] -= M*matrix1[k][j];
	}

	for (i = 0; i < n; i++)
	{
		det *= matrix1[i][i];
	}
	for (i = 0; i<n; i++)
	{
		for (j = 0; j < n; j++)
			cout << matrix1[i][j] << " ";
		cout << endl;
	}
	return det;
}

float * sideDeters(float **matrix1, float * b, int n) {
	float  * deters = new float[n];
	int i, j;

	float **matrix = new float*[n];

	for (j = 0; j<n; j++)
	{
		mtrxCopy(matrix1, matrix, n, n);
		for (i = 0; i<n; i++)
		{
			matrix[i][j] = b[i];
		}
		deters[j] = determinant(matrix, n);
		
	}
	return deters;
}

float * methodGauss(float **matrix1, int n, int m) {
	int i, j;
	float **matrix = new float*[n];
	mtrxCopy(matrix1, matrix, n, m);

	//Метод Гаусса
	//Прямой ход,приведение матрицы к треугольному виду
	float  str, *back = new float[m];
	int k;
	for (i = 0; i<n; i++)
	{
		str = matrix[i][i];
		for (j = m-1; j >= i; j--)
			matrix[i][j] /= str;
		for (j = i + 1; j<n; j++)
		{
			str = matrix[j][i];
			for (k = m-1; k >= i; k--)
				matrix[j][k] -= str*matrix[i][k];
		}
	}

	cout << "Матрица преведенная к треугольному виду: \t " << endl;
	display(matrix, n, m);
	cout << endl;
	//Обратный ход,исключения
	back[n - 1] = matrix[n - 1][n];
	for (i = n - 2; i >= 0; i--)
	{
		back[i] = matrix[i][n];
		for (j = i + 1; j < n; j++)
			back[i] -= matrix[i][j] * back[j];
	}
	return back;
}
//выведение матрицы
void display(float **matrix, int n, int m) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			cout << matrix[i][j] << "  ";
		cout << endl;
	}
}
//копирование элементов в новую структуру
void mtrxCopy(float ** matrix, float ** matrix1, int n, int m) {
	int i, j;
	for (i = 0; i < n; i++)
		matrix1[i] = new float[m];

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			matrix1[i][j] = matrix[i][j];
		}
	}
}

//транспонирование матрицы
void TransponMtx(float **matrix, float **tMatr, int n) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			tMatr[j][i] = matrix[i][j];
}
//освобождения памяти
void FreeMem(float **matrix, int n)
{
	for (int i = 0; i < n; i++)
		delete[] matrix[i];
	delete[] matrix;
}
//печати матрицы
void PrintMtx(float **matrix, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}
//вычеркивания строки и столбца
void Get_matr(float **matrix, int n, float **temp_matrix, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matrix[ki][kj] = matrix[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}
//вычисления определителя матрицы
int Det(float **matrix, int n)
{
	int determ = 0;   //временная переменная для хранения определителя
	int k = 1;      //степень
	if (n < 1) {
		cout << "Не верный размер матрицы!!!" << endl;
		return 0;
	}
	else if (n == 1)
		determ = matrix[0][0];
	else if (n == 2)
		determ = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
	else {
		for (int i = 0; i < n; i++) {
			int m = n - 1;
			float **temp_matr = new float *[m];
			for (int j = 0; j < m; j++)
				temp_matr[j] = new float[m];
			Get_matr(matrix, n, temp_matr, 0, i);
			determ = determ + k * matrix[0][i] * Det(temp_matr, m);
			k = -k;
			FreeMem(temp_matr, m);
		}
	}
	return determ;
}
