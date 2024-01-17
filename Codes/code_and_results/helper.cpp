#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <set>
#include <stdlib.h>
using namespace std;

int num_frames, num_atoms, num_coordinates, num_residues;
string boxes_data_file = "boxes.txt";
string pdb_file_name = "input.pdb";

class Box_Stats{
    public:
        double x_coord, y_coord, z_coord;
    
        void print_stats(){
            cout << "Coordinates: " << x_coord << ", " << y_coord << ", " << z_coord << '\n';
            cout << '\n';
        }
};

// Function to extract coordinates from the pdb file
class Molecule_Stats{
    public:
        int atom_index, molecule_index;
        double x_coord, y_coord, z_coord;
        string character_1, character_2;
    
        void print_stats(){
            cout << "Atom Index: " << atom_index << '\n';
            cout << "Molecule Index: " << molecule_index << '\n';
            cout << "Coordinates: " << x_coord << ", " << y_coord << ", " << z_coord << '\n';
            cout << "Atom Name: " << character_1 << '\n';
            cout << "Place : " << character_2 << '\n';
            cout << '\n';
        }
};

int get_residues(vector<Molecule_Stats> &molecule_list){
    set<int> residues;
    for (auto &item : molecule_list){
        residues.insert(item.molecule_index);
    }
    return residues.size();
}

vector<Box_Stats> box_coordinates;
vector<Molecule_Stats> molecule_list;

void extractor(string &file_name) {
    molecule_list.clear();

    ifstream file_in(file_name);
    char * dump_line = (char *) malloc(150*sizeof(char));
    string first_word;
    double dump_double;
    // int line_no = 1;
    if(file_in){
        while(!file_in.eof() && file_in){
            file_in >> first_word;
            if(first_word != "ATOM"){
                file_in.getline(dump_line, 149);
            }
            else if(first_word == "ENDMDL"){
                break;
            }
            else{
                Molecule_Stats molecule;
                file_in >> molecule.atom_index;
                file_in >> molecule.character_1 >> molecule.character_2;
                file_in >> molecule.molecule_index;
                file_in >> molecule.x_coord >> molecule.y_coord >> molecule.z_coord;
                file_in >> dump_double >> dump_double;
                molecule_list.push_back(molecule);
                molecule.print_stats();
            }
        }
    }
    return;
}

void fill_box_vec() {
    ifstream file_in(boxes_data_file);
    char * dump_line = (char *) malloc(150*sizeof(char));
    int line_no = 1;
    if(file_in){
        while(!file_in.eof() && file_in){
            Box_Stats box;
            cin >> box.x_coord >> box.y_coord >> box.z_coord;
            box_coordinates.push_back(box);
        }
    }
    return;
}




int main() {
    num_frames = 1;
    num_coordinates = 3;
    extractor(pdb_file_name);
    num_atoms = molecule_list.size();
    num_residues = get_residues(molecule_list);
    // cout << molecule_list.size() << " :: Size\n";
    // for(auto &item: molecule_list){
    //     item.print_stats();
    // }
    system("python3 get_boxes.py");
    for(auto &item : box_coordinates){
        item.print_stats();
    }
    return 0;
}