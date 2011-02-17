#include <iostream>
#include <fstream>
#include <assert.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

using namespace std;
using namespace boost;

class Base :
    public boost::enable_shared_from_this<Base> ,
    private boost::noncopyable {
protected:
    Base() {}
    virtual ~Base() {}
public:
    virtual const std::string& name() const = 0;
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int){}
};
typedef boost::shared_ptr<Base> base_ptr;

class Derived : public Base {
public:
    Derived(const std::string& name )
               : Base(), name_(name){}
    Derived() {}
    virtual ~Derived() {}

    virtual const std::string& name() const { return name_;}

    base_ptr find(const std::string& name) const {
        if (name_ == name) {
            Derived* nonConstThis = const_cast<Derived*>(this);
            return nonConstThis->shared_from_this(); // throws after restore
        }
        return base_ptr();
    }

private:
    std::string name_;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int) {
        ar & boost::serialization::base_object<Base>(*this);
        ar & name_;
    }
};
typedef boost::shared_ptr<Derived> d_ptr;

class Holder {
public:
    Holder() {}
    void add(d_ptr b) { vec_.push_back(b);}

    base_ptr find(const std::string& name) const {
        for(size_t i=0; i < vec_.size(); i++){
            base_ptr base = vec_[i]->find(name);
            if (base) return base;
        }
        return base_ptr();
    }
private:
    std::vector<d_ptr> vec_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int){
         ar & vec_;
    }
};

void save(const Holder& holder, const char* filename)
{
    std::ofstream ofs( filename );
    boost::archive::text_oarchive oa( ofs );
    oa << holder;
}

void restore(Holder& holder, const char* filename)
{
    std::ifstream ifs( filename );
    boost::archive::text_iarchive ia( ifs );
    ia >> holder;
}

int main()
{
    std::string fileName = "test.txt";

    {
        Holder holder;
        holder.add( d_ptr(new Derived("me")));
        holder.add( d_ptr(new Derived("you")));

        base_ptr b = holder.find("me");
        assert(b.get()); // works

        save(holder, fileName.c_str());
    }

    Holder restored;
    restore(restored,fileName.c_str());

    try {
        base_ptr b = restored.find("me"); // fails ?
        assert(b.get()); // never gets here
    }
    catch (std::exception& e) {
        std::cout << "Ooops exception " << e.what()
        << " thrown, test failed??? \n";
        std::remove(fileName.c_str());
        return 1;
    }

    //std::remove(fileName.c_str());
    cout << "test passed\n";
} 
