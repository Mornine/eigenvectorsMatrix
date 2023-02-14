#include "holder.h"
#include "Givens.h"
using namespace std;
/* максимальное число итераций */
#define MAXITER 1000

/* Здесь определяются некоторые утилиты типа выделения памяти */
void tred2(double **a, int n, double *d, double *e) {
    int l, k, j, i;
    double scale, hh, h, g, f;

    for (i = n; i >= 2; i--) {
        l = i - 1;
        h = scale = 0.;
        if (l > 1) {
            /* вычислить шкалу */
            for (k = 1; k <= l; k++) scale += fabs(a[i][k]);
            /* малая величина шкалы -> пропустить преобразование */
            if (scale == 0.) e[i] = a[i][l];
            else {
                /* отмасштабировать строку и вычислить s2 в h */
                for (k = 1; k <= l; k++) {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }
                /* вычислить вектор u */
                f = a[i][l];
                g = (f >= 0. ? -sqrt(h) : sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                /* записать u на место i-го ряда a */
                a[i][l] = f - g;
                /* вычисление u/h, Au, p, K */
                f = 0.;
                for (j = 1; j <= l; j++) {
                    a[j][i] = a[i][j] / h;
                    /* сформировать элемент Au (в g) */
                    g = 0.;
                    for (k = 1; k <= j; k++) g += a[j][k] * a[i][k];
                    for (k = j + 1; k <= l; k++) g += a[k][j] * a[i][k];
                    /* загрузить элемент p во временно неиспользуемую область e */
                    e[j] = g / h;
                    /* подготовка к формированию K */
                    f += e[j] * a[i][j];
                }
                /* Сформировать K */
                hh = f / (h + h);
                for (j = 1; j <= l; j++) {
                    /* Сформировать q и поместить на место p (в e) */
                    f = a[i][j];
                    e[j] = g = e[j] - hh * f;
                    /* Трансформировать матрицу a */
                    for (k = 1; k <= j; k++) a[j][k] -= (f * e[k] + g * a[i][k]);
                }
            }
        } else e[i] = a[i][l];
        d[i] = h;
    }
    d[1] = 0.;
    e[1] = 0.;
    for (i = 1; i <= n; i++) {
        l = i - 1;

        if (d[i] != 0.) {
            for (j = 1; j <= l; j++) {
                g = 0.;
                /* формируем PQ, используя u и u/H */
                for (k = 1; k <= l; k++) g += a[i][k] * a[k][j];
                for (k = 1; k <= l; k++) a[k][j] -= g * a[k][i];
            }
        }
        d[i] = a[i][i];

        a[i][i] = 1;
        for (j = 1; j <= l; j++) a[j][i] = a[i][j] = 0.;
    }
}


/* максимальное число итераций */
#define MAXITER 1000

void tqli(double *d, double *e, int n, double **z) {
    int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;

    for (i = 2; i <= n; i++) e[i - 1] = e[i];
    e[n] = 0.;

    for (l = 1; l <= n; l++) {

        iter = 0;

        do {

            for (m = l; m <= n - 1; m++) {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if ((float) fabs(e[m] + dd) == dd) break;
            }

            if (m != l) {

                if (++iter >= MAXITER) {
                   // cout << "Error!!!";
                    //cout << "диагностика " << endl;
                    return;
                };

                g = (d[l + 1] - d[l]) / (2. * e[l]);
                r = hypot(1., g);
                /* здесь d_m - k_s */
                if (g >= 0.) g += fabs(r);
                else g -= fabs(r);
                g = d[m] - d[l] + e[l] / g;

                s = c = 1.;
                p = 0.;

                for (i = m - 1; i >= l; i--) {
                    f = s * e[i];
                    b = c * e[i];
                    e[i + 1] = r = hypot(f, g);

                    if (r == 0.) {
                        d[i + 1] -= p;
                        e[m] = 0.;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2. * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;
                    /* Содержимое следующего ниже цикла необходимо опустить, если
                       не требуются значения собственных векторов */
                    for (k = 1; k <= n; k++) {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }

                if (r == 0. && i >= l) continue;

                d[l] -= p;
                e[l] = g;
                e[m] = 0.;
            }
        } while (m != l);
    }
}

int startHolder() {
    int n;
    int i, j, size, m;
    double z;
    double **a, *e, *d;

    ifstream input("input.txt");
    ofstream fout("outputHolder.txt", ios::out);

    fout << "Поиск собственных чисел и векторов матриц методом отражений Хаусхолдера.\n";

    i = 0;  //количество чисел в файле
    double temp; //временная переменная

    while (!input.eof()) //cчитываем до конца файла
    {
        if (!(input >> temp)) {
            fout << "\nОшибка! Выход за пределы диапазона.Проверьте, что элементы введенной матрицы являются числами,попадающими в диапазон от -1.7*10^308 до 1.7*10^308.\n";
            error();
            return 0;
        }
        i++;
    }

    //поиск количества чисел в одной строке
    input.seekg(0, ios::beg);
    input.clear();
    j = 0;
    char symbol;
    while (!input.eof()) {
        input.get(symbol);
        if (symbol == ' ') j++;
        if (symbol == '\n') break;     //если конец строки, то выход их цикла
    }

    input.seekg(0, ios::beg);
    input.clear();

    n = i / (j + 1); //число строк
    m = j + 1; //число столбцов

    if (m != n) {
        fout << "\nОшибка! Считываемая матрица не является квадратной.\nПроверьте, чтобы число строк матрицы было равно числу столбцов \n";
        error();
        return 0;
    }

    a = new double *[n + 1];
    for (i = 1; i < n + 1; i++) a[i] = new double[n + 1];

    for (i = 1; i < n + 1; i++) {
        for (j = 1; j < n + 1; j++) {
            input >> z;
            a[i][j] = z;
        }
    }

    fout << "\nРазмерность матрицы: " << m << "\n";

    fout << "\nМатрица А: \n";
    for (i = 1; i < n + 1; i++) {
        fout << i + 1 << " row: ";
        for (j = 1; j < n + 1; j++) {
            fout << setw(7) << a[i][j];
        }
        fout << "\n";
    }
    fout << "\n\n";
   d = new double[n];
    e = new double[n];
    for (i = 0; i < n; i++) d[i] = e[i] = 0;

    tred2(a, n, d, e);

    tqli(d, e, n, a);


    for (i = 1; i < n + 1; i++) {
        fout << "\nСобственный вектор k" << i + 1 << ":\n";
        for (j = 1; j < n + 1; j++) {
            fout << round(a[i][j] * 1000000) / 1000000 << "\n";
        }
    }

    fout << "\nСобственные значения:\n";
    for (i = 1; i < n + 1; i++) {
        fout << i << ") " << round(d[i] * 1000000) / 1000000 << "\n";
    }

    fout.close();
    input.close();

    cout << " Программа выполнена успешно.\n Результаты выполнения программы  были записаны в файл output.txt\n"
         << endl;
    cout << "\n\n Для закрытия консоли дважды нажмите Enter...";
    getchar();

    //for (i = 0; i < n; i++) delete a[i]; // удаление динамического массива
    return 0;
}
