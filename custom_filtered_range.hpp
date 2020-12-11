#include <boost/range/adaptor/argument_fwd.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/transform_iterator.hpp>

template<typename Value>
class replace_value
{
public:
    typedef const Value& result_type;
    typedef const Value& argument_type;

    replace_value(const Value& from, const Value& to)
        : m_from(from), m_to(to)
    {
    }

    const Value& operator()(const Value& x) const
    {
        return (x == m_from) ? m_to : x;
    }
private:
    Value m_from;
    Value m_to;
};

template<typename Fn, typename Range>
struct custom_filtered_range :
    public boost::iterator_range<
        boost::filter_iterator<
            Fn,
            typename boost::range_iterator<Range>::type
        >
    >
{
private:
    typedef typename boost::range_value<Range>::type value_type;
    typedef typename boost::range_iterator<Range>::type iterator_base;
    //typedef replace_value<value_type> Fn;
    typedef boost::filter_iterator<Fn, iterator_base> custom_filtered_iterator;
    typedef boost::iterator_range<custom_filtered_iterator> base_t;

public:
    custom_filtered_range(Fn f, Range& rng)
        : base_t(make_filter_iterator( f, boost::begin(rng), boost::end(rng)),
                 make_filter_iterator( f, boost::end(rng), boost::end(rng)))
     { }
};

template<typename T>
struct custom_filter_holder : public boost::range_detail::holder<T>
{
public:
    custom_filter_holder(T r) : boost::range_detail::holder<T>(r)
    { }
};

static boost::range_detail::forwarder<custom_filter_holder>
custom_filtered = boost::range_detail::forwarder<custom_filter_holder>();

template<typename Fn, typename SinglePassRange>
inline custom_filtered_range<Fn, SinglePassRange>
operator|(SinglePassRange& rng,
          const custom_filter_holder<Fn>& f)
{
    return custom_filtered_range<Fn, SinglePassRange>(f.val, rng);
}

template<typename Fn, typename SinglePassRange>
inline custom_filtered_range<Fn, const SinglePassRange>
operator|(const SinglePassRange& rng,
          const custom_filter_holder<Fn>& f)
{
    return custom_filtered_range<Fn, const SinglePassRange>(f.val, rng);
}
