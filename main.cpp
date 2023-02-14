#include "Givens.h"
#include "holder.h"

int main() {
    setlocale(LC_ALL, "rus");
    thread thread1(startGivens);
    thread1.join();  //ожидание завершения потока
    thread thread2(startHolder);
    thread2.join();  //ожидание завершения потока
    return 0;
}
