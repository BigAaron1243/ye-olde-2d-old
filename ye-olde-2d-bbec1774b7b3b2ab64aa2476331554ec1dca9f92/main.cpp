#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>

class Wall {
    public:
        std::vector<double> a, b;
        void set_values(double, double, double, double, int);
        std::vector<std::vector<float>> edge;
};

void Wall::set_values (double ax, double ay, double bx, double by, int setresolution) {
    double length_xy[2] = {0,0};
    length_xy[0] = bx - ax;
    length_xy[1] = by - ay;
    for (int i = 0; i < setresolution; i++) {
        edge.push_back({((length_xy[0] / setresolution) * i), ((length_xy[1] / setresolution) * i)});
    }
    a = {ax, ay};
    b = {bx, by};
}

int rounder(float d)
{
  return floor(d + 0.5);
}

class Manifold {
    public:
        double angle, primary_circle, secondary_circle, v1, v2;
        std::vector<double> Resultant_vector_primary();
};

std::vector<double> Manifold::Resultant_vector_primary() {

}

class Circle {
    public:
        //radius, x-coordinate, y-coordinate, x-velocity (m/s), y-velocity (m/s), mass (kg), coefficient of friction (mu), perimeter (m)
        double r, x, y, m, u, perimeter;
        double vx = 0;
        double vy = 0;
        double instantaneous_collision_angle, instantaneous_collision_velocity, instantaneous_displacement;
        bool enable_gravity, _static;
        bool collision = false;
        int collision_cooldown = 0;
        void set_values(double, double, double, double, int, bool, bool);
        std::vector<std::vector<float>> edge;
        void calculate_new_position(double, double);
        double velocity();
};

double Circle::velocity() {
    return sqrt(pow(vx, 2) + pow(vy, 2));
}

void Circle::set_values (double setr, double setx, double sety, double setm, int setresolution, bool set_gravity, bool set_static) {
    r = setr; x = setx; y = sety; m = setm;
    enable_gravity = set_gravity;
    _static = set_static;
    for (int i = 0; i < setresolution; i++) {
        edge.push_back({sin(i*((2*M_PI)/setresolution)) * setr,cos(i*((2*M_PI)/setresolution)) * setr});
    }
    perimeter = 2 * M_PI * setr;
}

void Circle::calculate_new_position(double delta_time, double ppm) {
    x = x + (ppm * vx * (delta_time / 1000000));
    y = y + (ppm * vy * (delta_time / 1000000));
}

void print_window (int plane[100][100]) {
    std::string canvas;
    for(size_t y = 100; y > 0; y--) {
        char line[200];
        for(size_t x = 0; x < 100; x++) {
            if(plane[x][y] == 1) {
                line[x * 2] = '0';
                line[x * 2 - 1] = '0';
            } else {
                line[x * 2] = ' ';
                line[x * 2 - 1] = ' ';
            }
        }
        canvas.append(line);
        canvas.append("\n");
    }
    std::cout << canvas << std::endl;
}

void set_console_cursor(int x, int y, HANDLE console) {
    COORD ccoord;
    ccoord.X = x;
    ccoord.Y = y;
    SetConsoleCursorPosition(console, ccoord);
}

double get_collision_angle_radians(double sum_radius, double distance_xy[2]) {
    double angle_of_bounce;
    double distance = pow(distance_xy[0], 2) + pow(distance_xy[1], 2);
    if (distance_xy[0] <= 0) {
        if (distance_xy[1] >= 0) {
            angle_of_bounce =  (abs(abs(atan(distance_xy[1] / distance_xy[0])) - M_PI));
        } else {
            angle_of_bounce =  (abs(atan(distance_xy[1] / distance_xy[0])) + (1 * M_PI));
        }
    }
    if (distance_xy[0] > 0) {
        if (distance_xy[1] >= 0) {
            angle_of_bounce =  (atan(distance_xy[1] / distance_xy[0]));
        } else {
            angle_of_bounce =  ((0.5 * M_PI) - abs(atan(distance_xy[1] / distance_xy[0]))) + (1.5 * M_PI);
        }
    }
    return angle_of_bounce;
}

int main()
{
    double pixels_per_metre = 5;
    int total_time;
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    std::vector<Circle> circle_list;
    std::vector<Wall> wall_list;
    Circle ball_3;
    ball_3.set_values(5,60,5,1,300, false, false);
    ball_3.vy = 5;
    circle_list.push_back(ball_3);
    Circle ball;
    ball.set_values(5,0,100,1,300, true, false);
    ball.vx = 5;
    circle_list.push_back(ball);
    Circle ball_2;
    ball_2.set_values(5,100,120,1,300, true, false);
    ball_2.vx = -5;
    circle_list.push_back(ball_2);

    Circle ball_1;
    ball_1.set_values(10,50,30,1,300, false, true);
    circle_list.push_back(ball_1);
    int canvas[100][100] = {0};

    while (true) {
        auto start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++) {
                canvas[i][j] = 0;
            }
        }
        for (size_t i = 0; i < wall_list.size(); i++) {
            for (size_t j = 0; j < wall_list[i].edge.size(); j++) {
                int xpos = rounder(wall_list[i].edge[j][0] + wall_list[i].a[0]);
                int ypos = rounder(wall_list[i].edge[j][1] + wall_list[i].a[1]);
                if (xpos < 100 && ypos < 100 && xpos > 0 && ypos > 0) {
                    canvas[xpos][ypos] = 1;
                }
            }
        }
        for (size_t i = 0; i < circle_list.size(); i++) {
            for (size_t j = 0; j < circle_list[i].edge.size(); j++) {
                int xpos = rounder(circle_list[i].edge[j][0] + circle_list[i].x);
                int ypos = rounder(circle_list[i].edge[j][1] + circle_list[i].y);
                if (xpos < 100 && ypos < 100 && xpos > 0 && ypos > 0) {
                    canvas[xpos][ypos] = 1;
                }
            }
        }
        //int lolcount;
        bool collision = false;
        for (size_t i = 0; i < circle_list.size(); i++) {
            for (size_t j = 0; j < circle_list.size(); j++) {
                double angle_of_bounce = 0;
                double distance_xy[2] = {0,0};
                distance_xy[0] = circle_list[i].x - circle_list[j].x;
                distance_xy[1] = circle_list[i].y - circle_list[j].y;
                double sum_radius = pow(circle_list[i].r + circle_list[j].r, 2);
                double distance = pow(distance_xy[0], 2) + pow(distance_xy[1], 2);
                if ((j != i) && (sum_radius > distance) && collision == false && circle_list[i].collision_cooldown < 1) {
                    circle_list[i].collision_cooldown = 10;
                    circle_list[j].collision_cooldown = 10;
                    circle_list[i].collision = true;
                    circle_list[j].collision = true;
                    //circle_list[i].instantaneous_displacement = sum_radius - distance;
                    //circle_list[j].instantaneous_displacement = sum_radius - distance;
                    circle_list[i].instantaneous_collision_velocity =  sqrt(pow(circle_list[j].vx, 2) + pow(circle_list[j].vy, 2)) + sqrt(pow(circle_list[i].vx, 2) + pow(circle_list[i].vy, 2)) * .5;
                    circle_list[j].instantaneous_collision_velocity =  sqrt(pow(circle_list[j].vx, 2) + pow(circle_list[j].vy, 2)) + sqrt(pow(circle_list[i].vx, 2) + pow(circle_list[i].vy, 2)) * .5;
                    double collision_angle = get_collision_angle_radians(sum_radius, distance_xy);
                    circle_list[i].instantaneous_collision_angle = collision_angle;
                    if (collision_angle > 1) {
                        circle_list[j].instantaneous_collision_angle = collision_angle - (M_PI);
                    } else {
                        circle_list[j].instantaneous_collision_angle = collision_angle + (M_PI);
                    }
                }
            }
        }

        set_console_cursor(0, 1, hConsole);
        print_window(canvas);
        auto stop = std::chrono::high_resolution_clock::now();
        auto delta_time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        total_time += delta_time.count();
        for (size_t i = 0; i < circle_list.size(); i++) {
            if (circle_list[i].enable_gravity == true) {
                circle_list[i].vy -= 9.8 * delta_time.count() / 1000000;
            }
            if (circle_list[i]._static == false && circle_list[i].collision == true) {
                circle_list[i].x += cos(circle_list[i].instantaneous_collision_angle) * circle_list[i].instantaneous_displacement;
                circle_list[i].y += sin(circle_list[i].instantaneous_collision_angle) * circle_list[i].instantaneous_displacement;
                circle_list[i].vx = cos(circle_list[i].instantaneous_collision_angle) * circle_list[i].instantaneous_collision_velocity;
                circle_list[i].vy = sin(circle_list[i].instantaneous_collision_angle) * circle_list[i].instantaneous_collision_velocity;
                circle_list[i].collision = false;
            }
            circle_list[i].calculate_new_position(delta_time.count(), pixels_per_metre);
            circle_list[i].collision_cooldown--;
        }

        set_console_cursor(0, 0, hConsole);
        std::cout << total_time << "      ";

    }
    return 0;
}
