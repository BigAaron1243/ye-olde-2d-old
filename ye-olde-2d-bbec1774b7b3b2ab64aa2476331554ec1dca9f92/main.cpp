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
        double get_angle();
        double get_gradient();
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

double Wall::get_gradient () {
return (b[1] - a[1]) / (b[0] - a[0]);
}


double Wall::get_angle () {
    return atan2((b[1] - a[1]), (b[0] - a[0]));
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
        double instantaneous_collision_angle, instantaneous_collision_velocity, instantaneous_displacement, normal_angle, reflection_angle;
        bool enable_gravity, _static;
        bool collision = false;
        int collision_cooldown = 0;
        void set_values(double, double, double, double, int, bool, bool);
        std::vector<std::vector<float>> edge;
        void calculate_new_position(double, double);
        double velocity();
        double vector_angle();
        int neg = 1;
};

double Circle::vector_angle() {
    double temp_angle = atan2(vy, vx);
    if (temp_angle > 0) {
        return temp_angle;
    } else {
        return temp_angle + (2 * M_PI);
    }

}

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

    Wall floor;
    floor.set_values(0, 5, 100, 5, 100);
    wall_list.push_back(floor);

    Circle ball_2;
    ball_2.set_values(5,82,90,0.5,300, false, false);
    ball_2.vy = -3;
    circle_list.push_back(ball_2);
    Circle ball_3;
    ball_3.set_values(5,60,5,1,300, false, false);
    ball_3.vy = 5;
    circle_list.push_back(ball_3);
    Circle ball;
    ball.set_values(5,0,100,1,300, true, false);
    ball.vx = 5;
    circle_list.push_back(ball);

    Circle ball_1;
    ball_1.set_values(10,50,30,3,300, false, false);
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
                    //circle_list[j].collision_cooldown = 10;
                    circle_list[i].collision = true;
                    //circle_list[j].collision = true;
                    //circle_list[i].instantaneous_displacement = sum_radius - distance;
                    //circle_list[j].instantaneous_displacement = sum_radius - distance;
                    circle_list[i].instantaneous_collision_velocity =  sqrt(pow(circle_list[j].vx, 2) + pow(circle_list[j].vy, 2)) + sqrt(pow(circle_list[i].vx, 2) + pow(circle_list[i].vy, 2)) * .5;
                    //circle_list[j].instantaneous_collision_velocity =  sqrt(pow(circle_list[j].vx, 2) + pow(circle_list[j].vy, 2)) + sqrt(pow(circle_list[i].vx, 2) + pow(circle_list[i].vy, 2)) * .5;
                    double collision_angle = get_collision_angle_radians(sum_radius, distance_xy);
                    circle_list[i].instantaneous_collision_angle = collision_angle;
                    set_console_cursor(210, 5, hConsole);
                    double b = (circle_list[i].vector_angle() - M_PI) - circle_list[i].instantaneous_collision_angle;
                    circle_list[i].reflection_angle = circle_list[i].instantaneous_collision_angle - b;
                    std::cout << circle_list[i].instantaneous_collision_angle << ", " << circle_list[i].vector_angle() << ", " << circle_list[i].vx << ", " << circle_list[i].vy << ", " << circle_list[i].reflection_angle << ", " << b;
                    set_console_cursor(circle_list[i].x * 2, circle_list[i].y * 2, hConsole);
                    std::cout << "+";

                    circle_list[i].neg = 1;
                    if (circle_list[i].vx > 0.01 && circle_list[j].vx > 0.01 ) {
                        //circle_list[i].reflection_angle = circle_list[i].reflection_angle + M_PI;
                        circle_list[i].neg = -1;
                    } else if (circle_list[i].vy > 0.01 && circle_list[j].vy > 0.01 ) {
                        //circle_list[i].reflection_angle = circle_list[i].reflection_angle + M_PI;
                        circle_list[i].neg = -1;
                    }
                    if (circle_list[i].vx < -0.01 && circle_list[j].vx < -0.01 ) {
                        //circle_list[i].reflection_angle = circle_list[i].reflection_angle + M_PI;
                        circle_list[i].neg = -1;
                    } else if (circle_list[i].vy < -0.01 && circle_list[j].vy < -0.01 ) {
                        //circle_list[i].reflection_angle = circle_list[i].reflection_angle + M_PI;
                        circle_list[i].neg = -1;
                    }

                    //_sleep(500);
                    //if (collision_angle > 1) {
                    //    circle_list[j].instantaneous_collision_angle = circle_list[i].instantaneous_collision_angle - (M_PI);
                    //} else {
                    //    circle_list[j].instantaneous_collision_angle = circle_list[i].instantaneous_collision_angle - (M_PI);
                    //}
                }
            }
           // for (size_t j = 0; j < wall_list.size(); j++) {
               // int setrays = 20;
                //for (int k = 0; k < setrays; k++) {
                    //if (((sin(k*((2*M_PI)/setrays)) * circle_list[i].r + circle_list[i].x - wall_list[j].a[0]) * wall_list[j].get_gradient() >  cos(k*((2*M_PI)/setrays)) * circle_list[i].r + circle_list[i].y -1 && sin(k*((2*M_PI)/setrays)) * circle_list[i].r + circle_list[i].x - wall_list[j].a[0]) * wall_list[j].get_gradient() <  cos(k*((2*M_PI)/setrays)) * circle_list[i].r + circle_list[i].y +1) {
                    //    set_console_cursor(210, 5, hConsole);
                    //    std::cout << "";
                    //}
                    //double wally = (circle_list[i].x - wall_list[j].a[0]) * wall_list[j].get_gradient();
                    //if (wally + 1 + (sin(k*((2*M_PI)/setrays)) * circle_list[i].r) > circle_list[i].y + (cos(k*((2*M_PI)/setrays)) * circle_list[i].r) && wally - 1 + (sin(k*((2*M_PI)/setrays)) * circle_list[i].r) < circle_list[i].y + (cos(k*((2*M_PI)/setrays)) * circle_list[i].r)) {
                    //    double b = circle_list[i].vector_angle() - M_PI - (wall_list[j].get_angle() + (M_PI * 0.5));
                    //    circle_list[i].reflection_angle = circle_list[i].instantaneous_collision_angle - b;
                    //}
               // }
            //}
        }

        set_console_cursor(0, 1, hConsole);
        print_window(canvas);
        auto stop = std::chrono::high_resolution_clock::now();
        auto delta_time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        total_time += delta_time.count();
        for (size_t i = 0; i < circle_list.size(); i++) {
            double movement_mod = 5;
            double speed_decay = 0.003;
            if (GetAsyncKeyState(VK_UP) < 0) {
                circle_list[0].vy += movement_mod * delta_time.count() / 1000000;
            } if (GetAsyncKeyState(VK_DOWN) < 0) {
                circle_list[0].vy -= movement_mod * delta_time.count() / 1000000;
            } if (GetAsyncKeyState(VK_LEFT) < 0) {
                circle_list[0].vx -= movement_mod * delta_time.count() / 1000000;
            } if (GetAsyncKeyState(VK_RIGHT) < 0) {
                circle_list[0].vx += movement_mod * delta_time.count() / 1000000;
            } if (GetAsyncKeyState(VK_LCONTROL) < 0) {
                circle_list[0].enable_gravity = true;
            }
            if (circle_list[i].enable_gravity == true) {
                circle_list[i].vy -= 9.8 * delta_time.count() / 1000000;
            }
            if (circle_list[i]._static == false && circle_list[i].collision == true) {
                circle_list[i].x += cos(circle_list[i].reflection_angle) * circle_list[i].instantaneous_displacement + 0.01;
                circle_list[i].y += sin(circle_list[i].reflection_angle) * circle_list[i].instantaneous_displacement + 0.01;
                circle_list[i].vx = cos(circle_list[i].reflection_angle) * circle_list[i].neg * circle_list[i].instantaneous_collision_velocity * 0.5 * (1/circle_list[i].m);
                circle_list[i].vy = sin(circle_list[i].reflection_angle) * circle_list[i].neg * circle_list[i].instantaneous_collision_velocity * 0.5 * (1/circle_list[i].m);
                circle_list[i].collision = false;
            }
            circle_list[i].calculate_new_position(delta_time.count(), pixels_per_metre);
            if (circle_list[i].x > 100) {circle_list[i].vx = -circle_list[i].vx;}
            if (circle_list[i].x < 0) {circle_list[i].vx = -circle_list[i].vx;}
            if (circle_list[i].y > 100) {circle_list[i].vy = -circle_list[i].vy;}
            if (circle_list[i].y < 0) {circle_list[i].vy = -circle_list[i].vy;}
            //if (circle_list[i].vx > 0) {circle_list[i].vx += speed_decay * circle_list[i].vx;}
            //if (circle_list[i].vx < 0) {circle_list[i].vx -= speed_decay * circle_list[i].vx;}
            //if (circle_list[i].vy > 0) {circle_list[i].vy += speed_decay * circle_list[i].vy;}
            //if (circle_list[i].vy < 0) {circle_list[i].vy -= speed_decay * circle_list[i].vy;}
            circle_list[i].collision_cooldown--;
        }

        set_console_cursor(0, 0, hConsole);
        std::cout << total_time << "      ";

    }
    return 0;
}
