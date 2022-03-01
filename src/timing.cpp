#include <gbp/timing.hpp>


Timer::Timer() : m_beg(clock_t::now())
{
}

void Timer::reset()
{
	m_beg = clock_t::now();
}

double Timer::elapsed() const
{
	return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
}

std::string Timer::contruction_time(std::string format)
{
	const std::time_t t_c = clock_t::to_time_t(m_beg);
	std::stringstream str_time;
	if (format == "asc")
	{
		str_time <<  std::put_time(std::localtime(&t_c), "%c %Z");
	}
	else if (format == "dir")
	{
		str_time << std::put_time(std::localtime(&t_c), "%m%d_%H%M%S");
	}
	return str_time.str();
}

std::string Timer::current_time(std::string format)
{
	const std::time_t t_c = clock_t::to_time_t(clock_t::now());
	std::stringstream str_time;
	if (format == "asc")
	{
		str_time <<  std::put_time(std::localtime(&t_c), "%c %Z");
	}
	else if (format == "dir")
	{
		str_time << std::put_time(std::localtime(&t_c), "%m%d_%H%M%S");
	}
	return str_time.str();
}