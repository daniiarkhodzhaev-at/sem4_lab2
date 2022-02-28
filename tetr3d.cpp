#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>
#include <octree_wd.c>

using namespace std;

struct Point {
    double x;
    double y;
    double z;
};

double pdist(const Point &a, const Point &b) {
    return std::sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

std::ostream &operator<<(std::ostream &o, const Point &p) {
    return o << '(' << p.x << ", " << p.y << ", " << p.z << ')';
}

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Скорость
    double vx;
    double vy;
    double vz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), vx(0.0), vy(0.0), vz(0.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double vx, double vy, double vz) 
            : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz)
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau) {
        x += vx * tau * 0;
        y += vy * tau * 0;
        z += vz * tau * 0;
    }
};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            // Модельная скалярная величина распределена как-то вот так
            nodes[i] = CalcNode(pointX, pointY, pointZ, 0.0, 0.0, 0.0);
        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau) {
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].move(tau);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        char x;
        std::string _;
        double t, x1, x2, x3, v1, v2, v3, p;
        std::ifstream fd("data/data_" + std::to_string(snap_number) + ".csv");
        std::getline(fd, _);
        std::vector<std::pair<Point, Point>> data;
        double minx, maxx, miny, maxy, minz, maxz;
        maxx = 0.0;
        maxy = 0.0;
        maxz = 0.0;
        minx = INFINITY;
        miny = INFINITY;
        minz = INFINITY;
        while (fd >> t >> x >> x1 >> x >> x2 >> x >> x3 >> x >> v1 >> x >> v2 >> x >> v3 >> x >> p) {
            // std::cout << t << x << x1 << x << x2 << x << x3 << x << v1 << x << v2 << x << v3 << '\n';
            data.push_back({{x1, x2, x3}, {v1, v2, v3}});
            minx = std::min(minx, x1);
            maxx = std::max(maxx, x1);
            miny = std::min(miny, x2);
            maxy = std::max(maxy, x2);
            minz = std::min(minz, x3);
            maxz = std::max(maxz, x3);
        }
        std::cout << data.size() << '\n';
        fd.close();

        double origin[OCTREE_COORD_COUNT] = {0.5*(minx + maxx), 0.5*(miny + maxy), 0.5*(minz + maxz)};
        double size[OCTREE_COORD_COUNT] = {1.0*(-minx + maxx)+1.0, 1.0*(-miny + maxy)+1.0, 1.0*(-minz + maxz)+1.0};
        Octree *root = makeTree(origin, size);
        for (auto x : data) {
            double pos[OCTREE_COORD_COUNT] = {x.first.x, x.first.y, x.first.z};
            double dat[OCTREE_COORD_COUNT] = {x.second.x, x.second.y, x.second.z};
            insert(root, pos, dat);
        }

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            Point _vel, x0, xx;
            double dist;
            dist = INFINITY;
            x0 = {nodes[i].x, nodes[i].y, nodes[i].z};

            /*
            for (const auto &xx : data) {
                if (pdist(xx.first, x0) < dist) {
                    dist = pdist(xx.first, x0);
                    _vel = xx.second;
                }
            }
            */
            Octree *t = NULL;
            double _x0[OCTREE_COORD_COUNT] = {x0.x, x0.y, x0.z,};
            Retcode rc = lookup(root, _x0, &t);
            if (t == NULL) {
                printf("No vertices found for vertex: %f %f %f\n", nodes[i].x, nodes[i].y, nodes[i].z);
                abort();
            }
            /* this means that there are no neighbours,
               we just simply need to go up a few levels
            if (t->pointCount == 0) {
                printf("No children for lookup result of %f %f %f\n", nodes[i].x, nodes[i].y, nodes[i].z);
                abort();
            }
            */
            while (t->pointCount == 0 && t->parentNode != NULL) {
                t = t->parentNode;
                /* looking for neighbours in children */
                for (size_t j = 0; j < OCTREE_NODE_COUNT; ++j) {
                    if (t->childNodes[j] != NULL && t->childNodes[j]->pointCount > 0) {
                        /* FIXME: here we will find not the nearest node, but it will work for now */
                        t = t->childNodes[j];
                        break;
                    }
                }
            }
            /* now this means that there are really no neighbours
               and this situation doesn't make any sense, so the only thing to do --- abort */
            if (t->pointCount == 0) {
                printf("No children for lookup result of %f %f %f\n", nodes[i].x, nodes[i].y, nodes[i].z);
                abort();
            }
            for (size_t j = 0; j < t->pointCount; ++j) {
                xx = {t->points[OCTREE_COORD_COUNT * j + 0],
                      t->points[OCTREE_COORD_COUNT * j + 1],
                      t->points[OCTREE_COORD_COUNT * j + 2],};
                if (pdist(xx, x0) < dist) {
                    dist = pdist(xx, x0);
                    _vel = {t->vecData[OCTREE_COORD_COUNT * j + 0],
                            t->vecData[OCTREE_COORD_COUNT * j + 1],
                            t->vecData[OCTREE_COORD_COUNT * j + 2],};
                }
            }

            // Добавляем значение векторного поля в этой точке
            double __vec[3] = {_vel.x, _vel.y, _vel.z};
            vel->InsertNextTuple(__vec);
        }

        std::cerr << "Finished octree manipulations\n";

        remove_octree(root);

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "tetr3d-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    // Шаг точек по пространству
    double h = 4.0;
    // Шаг по времени
    double tau = 0.01;

    const unsigned int GMSH_TETR_CODE = 4;

    // Теперь придётся немного упороться:
    // (а) построением сетки средствами gmsh,
    // (б) извлечением данных этой сетки в свой код.
    gmsh::initialize();
    gmsh::model::add("t13");

    // Считаем STL
    try {
        gmsh::merge("suzanne_drag.msh");
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    /*
    // Восстановим геометрию
    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    // gmsh::model::mesh::createGeometry();

    // Зададим объём по считанной поверхности
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    // Зададим мелкость желаемой сетки
    // int f = gmsh::model::mesh::field::add("MathEval");
    // gmsh::model::mesh::field::setString(f, "F", "4");
    // gmsh::model::mesh::field::setAsBackgroundMesh(f);

    // Построим сетку
    gmsh::model::mesh::generate(3);

    */
    gmsh::model::mesh::refine();
    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    // TODO: неплохо бы полноценно данные сетки проверять, да

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    size_t data_size = 60;
    for (size_t i = 0; i < data_size; ++i) {
        mesh.snapshot(i);
    }

    return 0;
}
