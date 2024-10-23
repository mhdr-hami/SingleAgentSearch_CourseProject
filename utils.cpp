#include "bits_stdc++.h"
using namespace std;

int **make_array(const string &stp) {
    stringstream ss(stp);

    int **arr = new int *[4];
    for (int i = 0; i < 4; i++)
        arr[i] = new int[4];

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            string num;
            ss >> num;
            arr[i][j] = stoi(num);
        }
    }

    return arr;
}


int **pick_next_stp(int &cnt) {
    cnt += 1;
    int ** stp;

    fstream new_file;

    new_file.open("/Users/mohammadrezahami/Documents/University/SingleAgent/SAS-Project/100Korf.txt",
                  ios::in); //open a file to perform read operation using file object

    if (new_file.is_open()) {
        string tp;
        for (int i = 0; i < cnt - 1; i++) getline(new_file, tp);
        getline(new_file, tp);
        string tp2;
        for (int i = 5; i < tp.size(); i++)
            tp2 += tp[i];

//        cout << tp2 << endl;
        return make_array(tp2);
    }
    return stp;
}

int *make_pancake(const string &stp) {
    stringstream ss(stp);

    int *arr = new int [16];

    for (int i = 0; i < 16; i++) {
            string num;
            ss >> num;
//            cout<<"num is "<<num<<endl;
            arr[i] = stoi(num);
    }

    return arr;
}

int *pick_next_pancake(int &cnt) {
    cnt += 1;
    int * stp;

    fstream new_file;

    new_file.open("/Users/mohammadrezahami/Documents/University/SingleAgent/SAS-Project/1000Pancakes_16.txt",
                  ios::in); //open a file to perform read operation using file object

    if (new_file.is_open()) {
        string tp;
        for (int i = 0; i < cnt - 1; i++) getline(new_file, tp);
        getline(new_file, tp);
        string tp2;
        for (int i = 0; i < tp.size(); i++)
            tp2 += tp[i];

//        cout << "tp2 is " <<tp2 << endl;
        return make_pancake(tp2);
    }
    return stp;
}

//4 3 5 6 2 1 0 7 8 9 10 11
//8 7 6 3 4 5 10 9 0 1 2 11
//4 3 6 7 8 5 10 9 0 1 2 11
//2 1 0 9 10 5 8 7 6 3 4 11
//11 4 3 6 7 8 5 10 9 0 1 2
//5 8 7 6 3 4 11 10 9 0 1 2
//7 8 5 6 3 4 11 10 9 0 1 2