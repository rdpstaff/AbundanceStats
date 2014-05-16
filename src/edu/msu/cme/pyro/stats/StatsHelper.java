/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.pyro.stats;

import java.util.List;

/**
 *
 * @author fishjord
 */
public class StatsHelper {
    public static <T extends Number> double getAverage(List<T> numbers) {
        double avg = 0;
        for(T i : numbers)
            avg += i.doubleValue();

        return avg / numbers.size();
    }

    public static <T extends Number> double getStdDeviation(List<T> numbers) {
        return getStdDeviation(numbers, getAverage(numbers));
    }

    public static <T extends Number> double getStdDeviation(List<T> numbers, double average) {
        double deviation = 0;
        for(T i : numbers)
            deviation += Math.pow(i.doubleValue() - average, 2);

        return Math.sqrt(deviation / numbers.size());
    }
}
