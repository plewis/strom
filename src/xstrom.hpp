#pragma once

class XStrom : public std::exception
    {
    public:
                            XStrom() throw() {}
                            XStrom(const std::string s) throw() : _msg() {_msg = s;}
        virtual             ~XStrom() throw() {}
        const char *        what() const throw() {return _msg.c_str();}
        
    private:
    
        std::string         _msg;
    }; 
