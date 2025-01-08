#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <functional>
#include <cfloat>
#include <climits>

using namespace std;

// Structure to represent a road segment
struct Road
{
    int start;
    int end;
    double capacity;   // Maximum vehicle capacity
    double traffic;    // Current traffic (number of vehicles)
    double base_rate;  // Base toll rate
    double toll_rate;  // Current toll rate
};

// Class for Congestion Graph
class CongestionGraph
{
private:
    int V;  // Number of places (nodes)
    vector<vector<pair<int, pair<int, float>>>> adjList;  // {destination, {weight, congestion_factor}}

public:
    CongestionGraph(int V);
    void addRoad(int u, int v, int w, float congestionFactor);
    void removeRoad(int u, int v);
    void displayNetwork();
    float findShortestPath(int src, int dest);
    float calculateTrafficCongestion();
    void updateRoadCongestion(int u, int v, float congestionFactor);
    vector<int> findAllPaths(int src, int dest);
    pair<int, float> getMaxCongestedRoad();
    bool isConnected();
    void displayCongestionInformation();
    int findBestPlaceToGo();
};

// Constructor
CongestionGraph::CongestionGraph(int V)
{
    this->V = V;
    adjList.resize(V);
}

// Add road with congestion factor
void CongestionGraph::addRoad(int u, int v, int w, float congestionFactor)
{
    adjList[u].push_back({v, {w, congestionFactor}});
    adjList[v].push_back({u, {w, congestionFactor}});
}

// Remove road between u and v
void CongestionGraph::removeRoad(int u, int v)
{
    adjList[u].erase(remove_if(adjList[u].begin(), adjList[u].end(),
                                [v](const pair<int, pair<int, float>>& edge) { return edge.first == v; }),
                      adjList[u].end());
    adjList[v].erase(remove_if(adjList[v].begin(), adjList[v].end(),
                                [u](const pair<int, pair<int, float>>& edge) { return edge.first == u; }),
                      adjList[v].end());
}

// Display the entire network
void CongestionGraph::displayNetwork()
{
    cout << "\nNetwork Connections:\n";
    for (int u = 0; u < V; u++) {
        cout << "\nPlace " << u << " is connected to:\n";
        for (const auto& edge : adjList[u])
        {
            int v = edge.first;
            int weight = edge.second.first;
            float congestionFactor = edge.second.second;
            cout << "  - Place " << v << " (Distance: " << weight << " km, Congestion: " << congestionFactor << ")\n";
        }
    }
}

// Find the shortest path considering congestion
float CongestionGraph::findShortestPath(int src, int dest)
{
    vector<float> dist(V, FLT_MAX);
    priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>> pq;
    pq.push({0.0, src});
    dist[src] = 0.0;

    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();

        for (const auto& neighbor : adjList[u])
        {
            int v = neighbor.first;
            int weight = neighbor.second.first;
            float congestionFactor = neighbor.second.second;
            float newWeight = weight * (1 + congestionFactor);

            if (dist[u] + newWeight < dist[v])
            {
                dist[v] = dist[u] + newWeight;
                pq.push({dist[v], v});
            }
        }
    }
    return dist[dest];
}

// Calculate the total traffic congestion in the network
float CongestionGraph::calculateTrafficCongestion()
{
    float totalCongestion = 0;
    for (int u = 0; u < V; u++)
    {
        for (const auto& edge : adjList[u])
        {
            totalCongestion += edge.second.second;
        }
    }
    return totalCongestion;
}

// Update the congestion factor for a road
void CongestionGraph::updateRoadCongestion(int u, int v, float congestionFactor)
{
    for (auto& edge : adjList[u])
    {
        if (edge.first == v)
        {
            edge.second.second = congestionFactor;
        }
    }
    for (auto& edge : adjList[v])
    {
        if (edge.first == u)
        {
            edge.second.second = congestionFactor;
        }
    }
}

// Find all paths between two nodes (using DFS or BFS)
vector<int> CongestionGraph::findAllPaths(int src, int dest)
{
    vector<int> paths;
    for (const auto& edge : adjList[src])
    {
        if (edge.first == dest)
        {
            paths.push_back(dest);
        }
    }
    return paths;
}

// Get the most congested road
pair<int, float> CongestionGraph::getMaxCongestedRoad()
{
    int maxRoad = -1;
    float maxCongestion = -1;
    for (int u = 0; u < V; u++)
    {
        for (const auto& edge : adjList[u])
        {
            if (edge.second.second > maxCongestion)
            {
                maxCongestion = edge.second.second;
                maxRoad = u;
            }
        }
    }
    return {maxRoad, maxCongestion};
}

// Check if the network is fully connected
bool CongestionGraph::isConnected()
{
    vector<bool> visited(V, false);
    queue<int> q;
    q.push(0);
    visited[0] = true;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (const auto& edge : adjList[u])
        {
            int v = edge.first;
            if (!visited[v])
            {
                visited[v] = true;
                q.push(v);
            }
        }
    }

    for (bool v : visited)
    {
        if (!v) return false;
    }
    return true;
}

// Display congestion information for each road
void CongestionGraph::displayCongestionInformation()
{
    cout << "\nCurrent Congestion Information for all roads:\n";
    for (int u = 0; u < V; u++) {
        for (const auto& edge : adjList[u])
        {
            int v = edge.first;
            float congestion = edge.second.second;
            cout << "Road from Place " << u << " to Place " << v << " has congestion factor: " << congestion << "\n";
        }
    }
}

// Find the best place to go with the least congestion
int CongestionGraph::findBestPlaceToGo()
{
    float minCongestion = FLT_MAX;
    int bestPlace = -1;

    for (int u = 0; u < V; u++)
    {
        float totalCongestionForPlace = 0;
        for (const auto& edge : adjList[u])
        {
            totalCongestionForPlace += edge.second.second;
        }

        if (totalCongestionForPlace < minCongestion)
        {
            minCongestion = totalCongestionForPlace;
            bestPlace = u;
        }
    }
    return bestPlace;
}

// Class for Traffic Management System
class TrafficManagementSystem
{
private:
    int rows, cols;
    vector<vector<int>> cityGrid;
    map<pair<int, int>, vector<pair<int, int>>> adjacencyList;

public:
    TrafficManagementSystem(vector<vector<int>> grid);
    void buildGraph();
    int bfsShortestPath(pair<int, int> start, pair<int, int> end);
    vector<pair<int, int>> detectCongestionPoints();
    void dfs(pair<int, int> start, set<pair<int, int>>& visited);
    vector<tuple<int, pair<pair<int, int>, pair<int, int>>>> kruskalMST();
    vector<int> dijkstra(const vector<vector<pair<int, int>>>& graph, int source);
};

// Constructor
TrafficManagementSystem::TrafficManagementSystem(vector<vector<int>> grid) : cityGrid(grid)
{
    rows = cityGrid.size();
    cols = cityGrid[0].size();
    buildGraph();
}

// Build graph from city grid
void TrafficManagementSystem::buildGraph()
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (cityGrid[i][j] == 1)
            {
                vector<pair<int, int>> neighbors;
                for (auto [dx, dy] : vector<pair<int, int>>{{-1, 0}, {1, 0}, {0, -1}, {0, 1}})
                {
                    int x = i + dx, y = j + dy;
                    if (x >= 0 && x < rows && y >= 0 && y < cols && cityGrid[x][y] == 1)
                    {
                        neighbors.push_back({x, y });
                    }
                }
                adjacencyList[{i, j}] = neighbors;
            }
        }
    }
}

// BFS to find the shortest path in an unweighted graph
int TrafficManagementSystem::bfsShortestPath(pair<int, int> start, pair<int, int> end)
{
    map<pair<int, int>, int> distances;
    for (auto &[node, _] : adjacencyList)
    {
        distances[node] = INT_MAX;
    }
    distances[start] = 0;

    queue<pair<int, int>> q;
    q.push(start);

    while (!q.empty())
    {
        auto current = q.front();
        q.pop();

        for (auto neighbor : adjacencyList[current])
        {
            if (distances[neighbor] == INT_MAX)
            {
                distances[neighbor] = distances[current] + 1;
                q.push(neighbor);
            }
        }
    }

    return (distances[end] == INT_MAX) ? -1 : distances[end];
}

// Detect intersections with potential traffic congestion
vector<pair<int, int>> TrafficManagementSystem::detectCongestionPoints()
{
    vector<pair<int, int>> congestionPoints;
    for (auto &[node, neighbors] : adjacencyList)
    {
        if (neighbors.size() > 3) // Arbitrary threshold for congestion
        {
            congestionPoints.push_back(node);
        }
    }
    return congestionPoints;
}

// Depth-First Search to explore the graph
void TrafficManagementSystem::dfs(pair<int, int> start, set<pair<int, int>>& visited)
{
    visited.insert(start);
    for (auto neighbor : adjacencyList[start])
    {
        if (visited.find(neighbor) == visited.end())
        {
            dfs(neighbor, visited);
        }
    }
}

// Kruskal's algorithm for Minimum Spanning Tree
vector<tuple<int, pair<pair<int, int>, pair<int, int>>>> TrafficManagementSystem::kruskalMST()
{
    vector<tuple<int, pair<pair<int, int>, pair<int, int>>>> edges;
    set<pair<pair<int, int>, pair<int, int>>> visitedEdges;

    for (auto &[node, neighbors] : adjacencyList)
    {
        for (auto neighbor : neighbors)
        {
            if (visitedEdges.find({node, neighbor}) == visitedEdges.end() &&
                visitedEdges.find({neighbor, node}) == visitedEdges.end())
            {
                edges.push_back(make_tuple(1, make_pair(node, neighbor)));
                visitedEdges.insert({node, neighbor});
            }
        }
    }

    sort(edges.begin(), edges.end());

    map<pair<int, int>, pair<int, int>> parent;
    for (auto &[node, _] : adjacencyList)
    {
        parent[node] = node;
    }

    function<pair<int, int>(pair<int, int>)> find = [&](pair<int, int> node)
    {
        if (parent[node] != node)
        {
            parent[node] = find(parent[node]);
        }
        return parent[node];
    };

    auto unionNodes = [&](pair<int, int> node1, pair<int, int> node2)
    {
        auto root1 = find(node1);
        auto root2 = find(node2);
        if (root1 != root2)
        {
            parent[root2] = root1;
        }
    };

    vector<tuple<int, pair<pair<int, int>, pair<int, int>>>> mst;
    for (auto &[weight, edge] : edges)
    {
        auto [node1, node2] = edge;
        if (find(node1) != find(node2))
        {
            unionNodes(node1, node2);
            mst.push_back({weight, edge});
        }
    }

    return mst;
}

// Dijkstra's Algorithm: Find shortest paths from source to all nodes
vector<int> TrafficManagementSystem::dijkstra(const vector<vector<pair<int, int>>>& graph, int source)
{
    int n = graph.size();
    vector<int> dist(n, INT_MAX);
    dist[source] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> Q;
    Q.push({0, source});

    while (!Q.empty())
    {
        int u = Q.top().second;
        Q.pop();

        for (auto &neighbor : graph[u])
        {
            int v = neighbor.first;
            int weight = neighbor.second;
            int alt = dist[u] + weight;

            if (alt < dist[v])
            {
                dist[v] = alt;
                Q.push({alt, v});
            }
        }
    }

    return dist;
}

// Smart Parking System class
class SmartParkingSystem
{
private:
    struct Car
    {
        string plateNumber;
        string carType; // For example, "Small", "Medium", "Large"
        int arrivalTime; // This could be a timestamp or just the time of day
    };

    struct ParkingSpot
    {
        int spotNumber;
        string spotSize; // "Small", "Medium", "Large"
        bool isOccupied;

        ParkingSpot(int number, string size) : spotNumber(number), spotSize(size), isOccupied(false) {}
    };

    struct CompareCar
    {
        bool operator()(const Car& a, const Car& b)
        {
            return a.arrivalTime > b.arrivalTime; // Prioritize earlier arriving cars
        }
    };

    vector<ParkingSpot> parkingLot;
    unordered_map<string, int> carInSpot; // Maps plate number to the spot number where the car is parked
    priority_queue<Car, vector<Car>, CompareCar> carQueue; // To simulate car arrivals

public:
    SmartParkingSystem(int totalSpots);
    void displayAvailableSpots();
    void displayParkingLotStatus();
    bool isValidCarType(char carType);
    string getCarTypeName(char carType);
    bool isValidPlateNumber(const string& plateNumber);
    void carArrives(const string& plateNumber, char carType, int arrivalTime);
    void assignParking(const Car& car);
    void carDeparts(const string& plateNumber);
    void showMenu();
    void interactiveMode();
};

// Constructor to initialize parking spots
SmartParkingSystem::SmartParkingSystem(int totalSpots)
{
    for (int i = 0; i < totalSpots; ++i)
    {
        string size = (i % 3 == 0) ? "Small" : (i % 3 == 1) ? "Medium" : "Large";
        parkingLot.push_back(ParkingSpot(i + 1, size));
    }
}

// Display available parking spots
void SmartParkingSystem::displayAvailableSpots()
{
    cout << "Available Parking Spots:\n";
    bool available = false;
    for (const auto& spot : parkingLot)
    {
        if (!spot.isOccupied)
        {
            cout << "Spot #" << spot.spotNumber << " (" << spot.spotSize << ")\n";
            available = true;
        }
    }
    if (!available)
    {
        cout << "No available spots.\n";
    }
}

// Show parking lot status after every action
void SmartParkingSystem::displayParkingLotStatus()
{
    cout << "Current Parking Lot Status:\n";
    for (const auto& spot : parkingLot)
    {
        string status = spot.isOccupied ? "Occupied" : "Available";
        cout << "Spot #" << spot.spotNumber << " (" << spot.spotSize << ") - " << status << "\n";
    }
}

// Check for valid car type
bool SmartParkingSystem::isValidCarType(char carType)
{
    return carType == 's' || carType == 'm' || carType == 'l';
}

// Map the car type to the full name
string SmartParkingSystem::getCarTypeName(char carType)
{
    if (carType == 's') return "Small";
    if (carType == 'm') return "Medium";
    if (carType == 'l') return "Large";
    return ""; // Default empty string if invalid
}

// Check if the plate number is valid
bool SmartParkingSystem::isValidPlateNumber(const string& plateNumber)
{
    if (plateNumber.empty())
    {
        return false;
    }
    for (char ch : plateNumber)
    {
        if (!isalnum(ch))
        {
            return false;
        }
    }
    return true;
}

// Car arrives and requests parking
void SmartParkingSystem::carArrives(const string& plateNumber, char carType, int arrivalTime)
{
    if (!isValidCarType(carType))
    {
        cout << "Invalid car type. Please enter 's' for Small, 'm' for Medium, or 'l' for Large.\n";
        return;
    }
    if (!isValidPlateNumber(plateNumber))
    {
        cout << "Invalid plate number. Plate number must be alphanumeric and cannot be empty.\n";
        return;
    }

    Car newCar = { plateNumber, getCarTypeName(carType), arrivalTime };
    carQueue.push(newCar);
    cout << "Car " << plateNumber << " arriving at time " << arrivalTime << "\n";
    assignParking(newCar);
}

// Assign parking to the arriving car
void SmartParkingSystem::assignParking(const Car& car)
{
    bool assigned = false;
    for (auto& spot : parkingLot)
    {
        if (!spot.isOccupied && (spot.spotSize == car.carType || spot.spotSize == "Large"))
        {
            spot.isOccupied = true;
            carInSpot[car.plateNumber] = spot.spotNumber;
            cout << "Car " << car.plateNumber
 << " parked in Spot #" << spot.spotNumber << " (" << spot.spotSize << ")\n";
            assigned = true;
            break;
        }
    }
    if (!assigned)
    {
        cout << "No available spot for car " << car.plateNumber << " (" << car.carType << "). All spots are full.\n";
    }
}

// Car departs and frees up the parking spot
void SmartParkingSystem::carDeparts(const string& plateNumber)
{
    auto it = carInSpot.find(plateNumber);
    if (it != carInSpot.end())
    {
        int spotNumber = it->second;
        parkingLot[spotNumber - 1].isOccupied = false;
        carInSpot.erase(it);
        cout << "Car " << plateNumber << " has left the parking. Spot #" << spotNumber << " is now available.\n";
    }
    else
    {
        cout << "No car found with plate number " << plateNumber << " in the parking lot.\n";
    }
}

// Show menu options for the user
void SmartParkingSystem::showMenu()
{
    cout << "1. Add Car\n";
    cout << "2. Remove Car\n";
    cout << "3. View Available Spots\n";
    cout << "4. View Parking Lot Status\n";
    cout << "5. Exit\n";
}

// Interactive mode for user to continuously interact with the system
void SmartParkingSystem::interactiveMode()
{
    int choice;
    string plateNumber;
    char carType;
    int arrivalTime;

    while (true)
    {
        showMenu();
        cout << "Enter your choice: ";
        cin >> choice;

        cin.ignore();  // Ignore the leftover newline character from previous input

        switch (choice)
        {
            case 1: // Add car
                cout << "Enter car plate number: ";
                getline(cin, plateNumber);
                cout << "Enter car type ('s' for Small, 'm' for Medium, 'l' for Large): ";
                cin >> carType;
                cout << "Enter arrival time (in minutes): ";
                cin >> arrivalTime;
                carArrives(plateNumber, carType, arrivalTime);
                break;
            case 2: // Remove car
                cout << "Enter car plate number to remove: ";
                cin.ignore();
                getline(cin, plateNumber);
                carDeparts(plateNumber);
                break;
            case 3: // View available spots
                displayAvailableSpots();
                break;
            case 4: // View parking lot status
                displayParkingLotStatus();
                break;
            case 5: // Exit
                cout << "Exiting system...\n";
                return;
            default:
                cout << "Invalid choice. Try again.\n";
                break;
        }
    }
}

// Function to adjust toll rates based on congestion
void adjustTollRates(vector<Road> &roads, double penalty_factor, double incentive_factor, double threshold)
{
    for (auto &road : roads)
    {
        double congestion_level = road.traffic / road.capacity;
        if (congestion_level > threshold)
        {
            road.toll_rate = road.base_rate + penalty_factor * (congestion_level - threshold);
        }
        else
        {
            road.toll_rate = max(road.base_rate - incentive_factor * (threshold - congestion_level), 0.0);
        }
    }
}

// Function to display the current toll rates
void displayTollRates(const vector<Road> &roads)
{
    cout << "\nUpdated Toll Rates for Each Road:\n";
    for (const auto &road : roads)
    {
        cout << "From Location " << road.start << " to Location " << road.end
             << ": Toll Rate = " << road.toll_rate << " currency units\n";
    }
}

// Helper function to get positive double input
double getPositiveDouble(const string &prompt)
{
    double value;
    while (true)
    {
        cout << prompt;
        cin >> value;
        if (cin.fail() || value < 0)
        {
            cout << "Invalid input. Please enter a non-negative number.\n";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        else
        {
            return value;
        }
    }
}

// Helper function to get double input within range
double getDoubleInRange(const string &prompt, double min, double max)
{
    double value;
    while (true)
    {
        cout << prompt;
        cin >> value;
        if (cin.fail() || value < min || value > max)
        {
            cout << "Invalid input. Please enter a number between " << min << " and " << max << ".\n";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\ n');
        }
        else
        {
            return value;
        }
    }
}

// Main function to run the entire system
int main()
{
    // Initialize Congestion Graph
    int V;
    cout << "Welcome to the Traffic Management System!\n";
    cout << "First, let's set up the number of places (nodes) in the city.\n";
    cout << "Enter the total number of places in your city: ";
    cin >> V;

    CongestionGraph cg(V);

    int roadCount;
    cout << "\nNow, let's add some roads! How many roads do you want to input? ";
    cin >> roadCount;

    cout << "\nPlease enter the details for each road:\n";
    for (int i = 0; i < roadCount; i++)
    {
        int u, v, w;
        float congestionFactor;
        cout << "\nEnter details for road " << i + 1 << " (Place1 Place2 Distance CongestionFactor): ";
        cin >> u >> v >> w >> congestionFactor;
        cg.addRoad(u, v, w, congestionFactor);
    }

    cg.displayNetwork();  // Display the network with congestion

    int src, dest;
    cout << "\nGreat! Now, let's calculate the shortest path between two places.\n";
    cout << "Enter the starting place and the destination place: ";
    cin >> src >> dest;
    float shortestPath = cg.findShortestPath(src, dest);
    cout << "The shortest path from Place " << src << " to Place " << dest << " considering traffic is: " << shortestPath << " km\n";

    cout << "\nTotal congestion in the city is: " << cg.calculateTrafficCongestion() << endl;

    int updateU, updateV;
    float newCongestion;
    cout << "\nWould you like to update the congestion factor of a road? (Enter '0' for no, '1' for yes): ";
    int updateChoice;
    cin >> updateChoice;
    if (updateChoice == 1)
    {
        cout << "Enter the road to update (Place1 Place2 NewCongestionFactor): ";
        cin >> updateU >> updateV >> newCongestion;
        cg.updateRoadCongestion(updateU, updateV, newCongestion);
    }

    cg.displayCongestionInformation();

    // Suggest the best place with minimum congestion
    int bestPlace = cg.findBestPlaceToGo();
    cout << "\nThe best place to go with the least congestion is Place " << bestPlace << ".\n";

    // Initialize Traffic Management System
    vector<vector<int>> cityGrid =
    {
        {1, 1, 0, 1},
        {1, 0, 1, 1},
        {1, 1, 1, 0},
        {0, 1, 0, 1}
    };

    TrafficManagementSystem trafficSystem(cityGrid);

    // Find shortest path
    auto shortestPathTraffic = trafficSystem.bfsShortestPath({0, 0}, {3, 3});
    cout << "Shortest Path: " << shortestPathTraffic << endl;

    // Detect congestion points
    auto congestionPoints = trafficSystem.detectCongestionPoints();
    cout << "Congestion Points:" << endl;
    for (auto point : congestionPoints)
    {
        cout << "(" << point.first << ", " << point.second << ")\n";
    }

    // DFS to explore the graph
    set<pair<int, int>> visited;
    trafficSystem.dfs({0, 0}, visited);
    cout << "Visited Nodes in DFS:" << endl;
    for (auto node : visited)
    {
        cout << "(" << node.first << ", " << node.second << ")\n";
    }

    // Find MST
    auto mst = trafficSystem.kruskalMST();
    cout << "Minimum Spanning Tree:" << endl;
    for (auto [weight, edge] : mst)
    {
        cout << "(" << get<0>(edge).first << ", " << get<0>(edge).second << ") -> (";
        cout << get<1>(edge).first << ", " << get<1>(edge).second << ")\n";
    }

    // Initialize Smart Parking System
    int totalSpots;
    cout << "Enter the total number of parking spots: ";
    cin >> totalSpots;

    SmartParkingSystem parkingSystem(totalSpots);
    parkingSystem.interactiveMode();

    // Initialize Roads for Toll Rate Adjustment
    int num_roads;
    cout << "Enter the number of road segments to monitor: ";
    cin >> num_roads;

    vector<Road> roads(num_roads);
    for (int i = 0; i < num_roads; ++i)
    {
        cout << "Enter details for Road " << i + 1 << "\n";
        cout << "  Starting location (numeric ID): ";
        cin >> roads[i].start;
        cout << "  Ending location (numeric ID): ";
        cin >> roads[i].end;
        roads[i].capacity = getPositiveDouble("  Maximum vehicle capacity for this road: ");
        roads[i].base_rate = getPositiveDouble("  Base toll rate (currency units): ");
        roads[i].traffic = 0; // Initialize traffic to 0
        roads[i].toll_rate = roads[i].base_rate;
    }

    // Input adjustment parameters
    double penalty_factor = getPositiveDouble("  Penalty factor for congestion (higher value increases toll rates faster): ");
    double incentive_factor = getPositiveDouble("  Incentive factor for free flow (higher value reduces toll rates faster): ");
    double congestion_threshold = getDoubleInRange("  Congestion threshold (as a fraction, e.g., 0.8 for 80% of capacity): ", 0.0, 1.0);

    char choice;
    do
    {
        // Update traffic data
        cout << "\nPlease enter the current traffic (number of vehicles) for each road:\n";
        for (int i = 0; i < num_roads; ++i)
        {
            roads[i].traffic = getPositiveDouble("  From Location " + to_string(roads[i].start) + " to Location " + to_string(roads[i].end) + ": ");
        }

        // Adjust toll rates
        adjustTollRates(roads, penalty_factor, incentive_factor, congestion_threshold);

        // Display updated toll rates
        displayTollRates(roads);

        // Ask if the user wants to continue
        cout << "\nDo you want to update traffic data again? Enter 'y' for yes or 'n' for no: ";
        cin >> choice;
        while (choice != 'y' && choice != 'Y' && choice != 'n' && choice != 'N')
        {
            cout << "Invalid input. Please enter 'y' for yes or 'n' for no: ";
            cin >> choice;
        }
    }
    while (choice == 'y' || choice == 'Y');

    cout << "\nThank you for using the Real-Time Toll Rate Adjustment System. Goodbye!\n";
    return 0;
}

/*
SAMPLE INPUT//////////////////////////////////////////////////////////////////////////////////////
Welcome to the Traffic Management System!
First, let's set up the number of places (nodes) in the city.
Enter the total number of places in your city: 4

Now, let's add some roads! How many roads do you want to input? 3

Please enter the details for each road:
Enter details for road 1 (Place1 Place2 Distance CongestionFactor): 0 1 10 0.2
Enter details for road 2 (Place1 Place2 Distance CongestionFactor): 1 2 5 0.1
Enter details for road 3 (Place1 Place2 Distance CongestionFactor): 2 3 15 0.3

Network Connections:

Place 0 is connected to:
  - Place 1 (Distance: 10 km, Congestion: 0.2)

Place 1 is connected to:
  - Place 0 (Distance: 10 km, Congestion: 0.2)
  - Place 2 (Distance: 5 km, Congestion: 0.1)

Place 2 is connected to:
  - Place 1 (Distance: 5 km, Congestion: 0.1)
  - Place 3 (Distance: 15 km, Congestion: 0.3)

Place 3 is connected to:
  - Place 2 (Distance: 15 km, Congestion: 0.3)

Great! Now, let's calculate the shortest path between two places.
Enter the starting place and the destination place: 0 3
The shortest path from Place 0 to Place 3 considering traffic is: 30 km

Total congestion in the city is: 0.6

Would you like to update the congestion factor of a road? (Enter '0' for no, '1' for yes): 1
Enter the road to update (Place1 Place2 NewCongestionFactor): 1 2 0.5

Current Congestion Information for all roads:
Road from Place 0 to Place 1 has congestion factor: 0.2
Road from Place 1 to Place 0 has congestion factor: 0.2
Road from Place 1 to Place 2 has congestion factor: 0.5
Road from Place 2 to Place 1 has congestion factor: 0.5
Road from Place 2 to Place 3 has congestion factor: 0.3
Road from Place 3 to Place 2 has congestion factor: 0.3

The best place to go with the least congestion is Place 0.

Enter the total number of parking spots: 5


SAMPLE output////////////////////////////////////////////////////////////////////////////////////////
Available Parking Spots:
Spot #1 (Small)
Spot #2 (Medium)
Spot #3 (Large)
Spot #4 (Small)
Spot #5 (Medium)

Current Parking Lot Status:
Spot #1 (Small) - Available
Spot #2 (Medium) - Available
Spot #3 (Large) - Available
Spot #4 (Small) - Available
Spot #5 (Medium) - Available

Enter your choice: 1
Enter car plate number: ABC123
Enter car type ('s' for Small, 'm' for Medium, 'l' for Large): s
Enter arrival time (in minutes): 10
Car ABC123 arriving at time 10
Car ABC123 parked in Spot #1 (Small)

Current Parking Lot Status:
Spot #1 (Small) - Occupied by ABC123
Spot #2 (Medium) - Available
Spot #3 (Large) - Available
Spot #4 (Small) - Available
Spot #5 (Medium) - Available

Enter your choice: 2
Car ABC123 has departed from spot 1.
/////////////////////////////////////////////////////////////////////////////////////////////////////
Explanation of the Code
Congestion Graph:

The CongestionGraph class represents a graph where nodes are places and edges are roads connecting those places.
The user can add roads with specific distances and congestion factors.
The graph can calculate the shortest path between two places, taking into account the congestion factors.
It can also calculate total congestion in the network and find the best place to go with the least congestion.
Traffic Management System:

The TrafficManagementSystem class represents a grid-based city layout.
It can build a graph from a 2D grid, find the shortest path using BFS, and detect congestion points based on the number of connections.
Smart Parking System:

The SmartParkingSystem class manages parking spots for cars.
It allows users to add cars, remove cars, and view available parking spots and the current status of the parking lot.
The system checks for valid car types and plate numbers, assigns parking spots based on availability, and maintains a record of which cars are parked where.
Toll Rate Adjustment:

The code includes functions to adjust toll rates based on current traffic conditions.
It calculates toll rates dynamically based on congestion levels and displays updated rates to the user.
Summary of the Interaction
The user starts by setting up the traffic management system, defining places and roads.
The user can then calculate the shortest path between two places and view congestion information.
After that, the user can interact with the smart parking system, adding and removing cars while checking the status of parking spots.
The system also allows for real-time adjustments of toll rates based on traffic data, providing a comprehensive traffic management solution.
This simulation demonstrates how the various components of the system work together to manage traffic and parking efficiently, providing users with valuable information and services.
*/
