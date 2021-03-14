#ifndef HEAT_LINALG_H
#define HEAT_LINALG_H

#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>
#include <stdexcept>

template <typename T>
class Vector{
public:
    // Default constructor
    Vector() :
            length(0), data(nullptr){}

    // Copy constructor
    Vector(const Vector& rhs) :
            length(rhs.length){
        data = new T[length];

        for (int i=0; i < length; i++)
        {
            data[i] = rhs.data[i];
        }
    }

    // Move constructor
    Vector(Vector&& rhs) noexcept :
            length(rhs.length), data(rhs.data){
        rhs.length = 0;
        rhs.data = nullptr;
    }

    // Size constructor
    explicit Vector(int n) :
            length(n), data(new T[n]){}

    // Initialiser list constructor
    [[maybe_unused]] Vector(std::initializer_list<T> l) :
            Vector(l.size()){
        std::uninitialized_copy(l.begin(),l.end(), data);
    }

    // Destructor
    ~Vector(){
        delete[] data;
        data = nullptr;
        length = 0;
    }

    // Copy assignment
    Vector& operator=(const Vector& rhs){
        if (this == &rhs){
            return *this;
        }

        delete[] data;
        data = new T[rhs.length];
        length = rhs.length;
        for (int i=0; i<length; i++)
        {
            data[i] = rhs.data[i];
        }
        return *this;
    }

    // Move assignment
    Vector& operator=(Vector&& rhs) noexcept {
        delete[] data;
        length = rhs.length;
        data = rhs.data;
        rhs.data = nullptr;
        rhs.length = 0;
        return *this;
    }

    // index operator
    T& operator[](int index){
        if (index > length - 1) { throw std::runtime_error("Index out of bounds."); }
        return data[index];
    }

    // const index operator
    const T& operator[](int index) const{
        if (index > length - 1) { throw std::runtime_error("Index out of bounds."); }
        return data[index];
    }

    // + operator
    template<typename U>
    auto operator+(const Vector<U>& rhs) const{
        if (length != rhs.length) {
            throw std::runtime_error("Addition of Vectors with different length.");
        }

        Vector< typename std::common_type<T,U>::type > result(length);
        for (int i=0; i<length; i++)
        {
            result[i] = data[i] + rhs[i];
        }

        return result;
    }

    // - operator
    template<typename U>
    auto operator-(const Vector<U>& rhs) const{
        if (length != rhs.length) {
            throw std::runtime_error("Subtraction of Vectors with different length.");
        }

        Vector< typename std::common_type<T,U>::type  > result(length);
        for (int i=0; i<length; i++)
        {
            result[i] = data[i] - rhs[i];
        }

        return result;
    }

    // scalar multiplication
    template<typename U>
    auto operator*(U scalar) const{
        Vector< typename std::common_type<T,U>::type  > result(length);
        for (int i=0; i<length; i++)
        {
            result[i] = data[i] * scalar;
        }

        return result;
    }

    // size of vector
    [[nodiscard]] int size() const{
        return length;
    }

    static Vector<T> ones(int n){
        Vector<T> v(n);
        for (int i=0; i<n; i++){
            v.data[i] = 1;
        }
        return v;
    }

private:
    int length{};
    T* data;
};

template<typename T>
std::ostream& operator<<(std::ostream& o, const Vector<T>& V){
    for (int i=0; i<V.size(); i++){
        o << V[i] << " ";
    }
    o << std::endl;
    return o;
}

template<typename T, typename U>
auto operator*(U scalar, const Vector<T>& lhs){
    Vector< typename std::common_type<T,U>::type  > result(lhs.size());
    for (int i=0; i<lhs.size(); i++)
    {
        result[i] = lhs[i] * scalar;
    }

    return result;
}

template<typename T, typename U>
typename std::common_type<T,U>::type
dot(const Vector<T>& lhs,
    const Vector<U>& rhs)
{
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("Dot of Vectors with different length.");
    }
    typename std::common_type<T,U>::type res = 0;
    for (int i=0; i<lhs.size(); i++)
    {
        res += lhs[i] * rhs[i];
    }
    return res;
}

template <typename T>
class Sparse
{
public:
    const std::pair<int,int> shape;
    std::map< std::pair<int,int>, T > internal;

    [[maybe_unused]] Sparse(int nrows, int ncols) :
            shape({nrows, ncols})
    {}

    ~Sparse(){
        internal.clear();
    }

    T& operator[](const std::pair<int,int>& ij){
        if ((ij.first > shape.first - 1) or (ij.second > shape.second - 1)) {
            throw std::runtime_error("Index out of bounds."); }

        return internal[ij];
    }

    const T& operator()(const std::pair<int, int>& ij) const{
        if (internal.find(ij) != internal.end()){
            return internal.at(ij);
        }
        else{
            throw std::runtime_error("Key not present");
        }
    }

    template<typename U>
    friend std::ostream& operator<<(std::ostream& o, const Sparse<U>& M);
};

template<typename T>
std::ostream& operator<<(std::ostream& o, const Sparse<T>& M)
{
    for (int i=0; i<M.shape.first; i++){

        for (int j=0; j<M.shape.second; j++){

            if (M.internal.find({i, j}) != M.internal.end()){
                std::cout << M.internal.at({i, j}) << " ";
            }
            else{
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
    }
    return o;
}

template<typename T, typename U>
Vector<typename std::common_type<T,U>::type>
operator*(const Sparse<T>& lhs,
          const Vector<U>& rhs)
{
    if (lhs.shape.second != rhs.size())
    {
        throw std::runtime_error("Incompatible shapes between matrix and vector");
    }

    Vector<typename std::common_type<T,U>::type> res(rhs.size());
    for (int i=0; i<rhs.size(); i++)
    {
        res[i] = 0;
    }

    for (auto const& it : lhs.internal)
    {
        int i = it.first.first;
        int j = it.first.second;

        res[i] += it.second * rhs[j];

    }

    return res;
}


#endif //HEAT_LINALG_H
