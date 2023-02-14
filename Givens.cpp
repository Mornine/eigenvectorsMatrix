#include "Givens.h"
void error() {
    std::cout << " Ошибка при выполнение программы. \n Текст ошибки записан в файл output.txt\n";
    std::cout << "\n\n Для закрытия консоли дважды нажмите Enter...";
    getchar();
}

bool isSimmetrial(double **coefficients, int number) {
    bool result = true;
    int i, j;
    for (i = 0; i < number; i++) {
        for (j = i + 1; j < number; j++) {
            if (coefficients[i][j] != coefficients[j][i]) {
                result = false;
                break;
            }
        }
        if (!result) { break; }
    }
    return result;
}


int rotation(double **coefficients, int number,
             double **solution, double precision) {
    int result = 1;
    int i, j, k;
    int maxI, maxJ;
    double max, fi;
    double **givensmatrix;
    givensmatrix = new double *[number]; //Выделение памяти под матрицу поворота
    for (i = 0; i < number; i++) {
        givensmatrix[i] = new double[number];
    }
    double **temp;                //Выделение памяти под временную переменную
    temp = new double *[number];
    for (i = 0; i < number; i++) {
        temp[i] = new double[number];
    }
    double fault = 0.0;
    for (i = 0; i < number; i++) {
        for (j = i + 1; j < number; j++) {
            fault = fault + coefficients[i][j] * coefficients[i][j];
        }
    }
    fault = sqrt(2 * fault);
    while (fault > precision) {
        max = 0.0;
        for (i = 0; i < number; i++) {             //Поиск максимального значения матрицы
            for (j = i + 1; j < number; j++) {
                if (coefficients[i][j] > 0 && coefficients[i][j] > max) {
                    max = coefficients[i][j];
                    maxI = i;
                    maxJ = j;
                } else if (coefficients[i][j] < 0 && -coefficients[i][j] > max) {
                    max = -coefficients[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }
        for (i = 0; i < number; i++) {               //Присваивание матрице поворота значения
            for (j = 0; j < number; j++) {
                givensmatrix[i][j] = 0;
            }
            givensmatrix[i][i] = 1;
        }  //Формирование матрицы вращения при повороте на прямой угол
        if (coefficients[maxI][maxI] == coefficients[maxJ][maxJ]) {
            givensmatrix[maxI][maxI] = givensmatrix[maxJ][maxJ] =
            givensmatrix[maxJ][maxI] = sqrt(2.0) / 2.0;
            givensmatrix[maxI][maxJ] = -sqrt(2.0) / 2.0;
        } else {
            fi = 0.5 * atan((2.0 * coefficients[maxI][maxJ]) /                                 // матрица вращения
                            (coefficients[maxI][maxI] - coefficients[maxJ][maxJ]));
            givensmatrix[maxI][maxI] = cos(fi);
            givensmatrix[maxJ][maxJ] = cos(fi);
            givensmatrix[maxI][maxJ] = -sin(fi);
            givensmatrix[maxJ][maxI] = sin(fi);

        }

        for (i = 0; i < number; i++) {   //Обнуление значения временной переменной
            for (j = 0; j < number; j++) {
                temp[i][j] = 0.0;
            }
        }
        for (i = 0; i < number; i++) {    //Формирование матрицы с поворотом элементов
            for (j = 0; j < number; j++) {
                for (k = 0; k < number; k++) {
                    temp[i][j] = temp[i][j] + givensmatrix[k][i] * coefficients[k][j];
                }
            }
        }
        for (i = 0; i < number; i++) {   //Обнуление исходной матрицы
            for (j = 0; j < number; j++) {
                coefficients[i][j] = 0.0;
            }
        }
        for (i = 0; i < number; i++) { //Присваивание исходной матрице значений,
            for (j = 0; j < number; j++) {    //полученных после поворота матрицы
                for (k = 0; k < number; k++) {
                    coefficients[i][j] = coefficients[i][j] +
                                         temp[i][k] * givensmatrix[k][j];
                }
            }
        }
        fault = 0.0;
        for (i = 0; i < number; i++) {          //Расчёт точности, для проверки условий цикла
            for (j = i + 1; j < number; j++) {
                fault = fault + coefficients[i][j] * coefficients[i][j];
            }
        }
        fault = sqrt(2 * fault);               //Обнуление временной переменной
        for (i = 0; i < number; i++) {
            for (j = 0; j < number; j++) {
                temp[i][j] = 0.0;
            }
        }
        for (i = 0; i < number; i++) {       //Формирование матрицы векторов во временной

            for (j = 0; j < number; j++) {
                for (k = 0; k < number; k++) {
                    temp[i][j] = temp[i][j] + solution[i][k] * givensmatrix[k][j];
                }
            }
        }
        for (i = 0; i < number; i++) {
            for (j = 0; j < number; j++) {         //Присваивание значений матрице векторов
                solution[i][j] = temp[i][j];//  из временной переменной
            }
        }
        result++;
    }
    return result;
}


int startGivens() {
    int i, j, m;
    int size;
    double **coefficients, **solution, precision = 0.001;
    ifstream input;
    ofstream fout;
    input.open("input.txt"); //открываем файл с исходными данными
    fout.open("outputGivens.txt");//открываем файл для записи результатов
    if (!fout.is_open() || !input.is_open()) {
        cout << "ошибка при открытии файлов!" << endl;
        return 0;
    }

    fout << "Поиск собственных чисел и векторов матриц методом вращений.\n";

    i = 0;  //число чисел в файле
    double temp; //временная переменная

    while (!input.eof()) //cчитываем до конца файла
    {
        if (!(input >> temp)) {
            fout<< "\nОшибка!Выход за пределы диапазона!\n Проверьте, что элементы введенной матрицы являются числами,попадающими в диапазон от -1.7*10^308 до 1.7*10^308.\n";
            error();
            return 0;
        }
        i++;
    }

    //поиск количества чисел в одной строке
    input.seekg(0, ios::beg);
    input.clear();

    //изначальное число пробелов
    j = 0;
    char symbol;
    while (!input.eof()) {
        input.get(symbol);
        if (symbol == ' ') j++;
        if (symbol == '\n') break;     //если конец строки, то выход их цикла
    }

    input.seekg(0, ios::beg);
    input.clear();

    size = i / (j + 1); //число строк
    m = j + 1; //число столбцов

    if (m != size) {
        fout
                << "\nОшибка! Считываемая матрица не является квадратной.\nПроверьте, чтобы число строк матрицы было равно числу столбцов матрицы.\n";
        error();
        return 0;
    }

    coefficients = new double *[size];
    solution = new double *[size];

    for (i = 0; i < size; i++) {
        coefficients[i] = new double[size];
        solution[i] = new double[size];
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            solution[i][j] = 0;
        }
        solution[i][i] = 1;
    }

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            input >> coefficients[i][j];
        }
    }


    //cout << "Введите точность расчета: ";                             //Ввод точности
    //cin >> precision;

    if (!isSimmetrial(coefficients, size)) {                           //Проверка на симметричность матрицы
        fout << "\nОшибка! Матрица не симметричная.\nПроверьте, что элементы матрицы симметричны относительно её диагонали.\n";
        error();
        return 0;
    } else {

        fout << "\nРазмерность матрицы: " << m << "\n";

        fout << "\nМатрица А: \n";
        for (i = 0; i < size; i++) {
            fout << i + 1 << " row: ";
            for (j = 0; j < size; j++) {
                fout << setw(7) << coefficients[i][j];
            }
            fout << "\n";
        }
        fout << "\n\n";

        int steps = rotation(coefficients, size, solution, precision);
        //cout << "\nРешение:\n\n";
        for (i = 0; i < size; i++) {

            double max = 0;

            //cout << "Норма вектора k" << i + 1 << ":\n";
            for (j = 0; j < size; j++) {
                if (solution[j][i] > max) {
                    max = solution[j][i];
                }
            }
           // max = pow(max, -1);
            fout << "\nСобственный вектор k" << i + 1 << ":\n";

            for (j = 0; j < size; j++) {
                solution[j][i] = solution[j][i] ;
                fout << solution[j][i] << "\n";
            }
            fout << " ";
        }

        fout << "\nСобственные значения:\n";
        for (i = 0; i < size; i++) {
            fout << i + 1 << ") " << coefficients[i][i] << "\n";
        }
        cout << "Общее число шагов: " << steps;   // Количество итераций
    }

    fout.close();
    input.close();

    cout << " Программа выполнена успешно.\n Результаты выполнения программы были записаны в файл outputGivens.txt\n"
         << endl;
 //cout << "\n\n Для закрытия консоли дважды нажмите Enter...";
   // getchar();

    return 0;
}
