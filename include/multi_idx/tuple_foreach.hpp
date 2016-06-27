#pragma once

template<std::size_t t_idx = 0, typename t_f, typename... t_tp>
inline typename std::enable_if<t_idx == sizeof...(t_tp), void>::type
tuple_foreach(std::tuple<t_tp...>& t, t_f& f)
{}

template<std::size_t t_idx = 0, typename t_f, typename... t_tp>
inline typename std::enable_if<t_idx < sizeof...(t_tp), void>::type
tuple_foreach(std::tuple<t_tp...>& t, t_f& f)
{
    f(std::get<t_idx>(t), t_idx);
    tuple_foreach<t_idx + 1, t_f, t_tp...>(t, f);
}

template<std::size_t t_idx = 0, typename t_f, typename... t_tp>
inline typename std::enable_if<t_idx == sizeof...(t_tp), void>::type
tuple_foreach(const std::tuple<t_tp...>& t, t_f& f)
{}

template<std::size_t t_idx = 0, typename t_f, typename... t_tp>
inline typename std::enable_if<t_idx < sizeof...(t_tp), void>::type
tuple_foreach(const std::tuple<t_tp...>& t, t_f& f)
{
    f(std::get<t_idx>(t), t_idx);
    tuple_foreach<t_idx + 1, t_f, t_tp...>(t, f);
}
