//
// Created by Work_PC on 28.12.2016.
//

#include "CalculateTrendParameters.h"

#include "Constants.h"
#include "MathUtils.h"
#include <algorithm>
#include <numeric>
#include "LinearTrend.h"

void addToTheEndOfArrayAndShift(int *array, int size, int valueToAdd);
void addToTheEndOfArrayAndShift(float *array, int size, float valueToAdd);

/**
 * ÐžÐ±Ð½ÑƒÐ»ÑÐµÑ‚ Ð¼Ð°ÑÑÐ¸Ð²
 * @param array ÑƒÐºÐ°Ð·Ð°Ñ‚ÐµÐ»ÑŒ Ð½Ð° Ð¼Ð°ÑÑÐ¸Ð²
 * @param size Ñ€Ð°Ð·Ð¼ÐµÑ€
 */
void reset(int array[], int size);

static const float InitAlphaType() {
    if (OPTIC_ALGORITHM_CONFIGURATION::Alstom) {
        return 0.07f;
    } else {
        return 0.1f;
    }

}

static const float InitPartOfAllowed() {
    if (OPTIC_ALGORITHM_CONFIGURATION::Alstom) {
        return 0.05;
    } else {
        return 0.15;
    }
}

static TrendParameters InitStartTrendParameters_ForVelocity(float alpha) {
    TrendParameters velocityTrendParameters;
    velocityTrendParameters.alpha = alpha;
    velocityTrendParameters.alpha1 = LinearTrend::calculateAlpha1(alpha);

    return velocityTrendParameters;
}

const int NumberOfOpticCalculationsInSecond = 12;

const int ErrorZeroArraySize = NumberOfOpticCalculationsInSecond * 15;
const int ErrorOtherArraySize = NumberOfOpticCalculationsInSecond * 2;
const int ShiftOtherArraySize = NumberOfOpticCalculationsInSecond * 2;

const float AlphaTypeConst = InitAlphaType();
const float PartOfAllowed = InitPartOfAllowed();

///!!!! Ð?ÐÐ?Ð¦Ð?ÐÐ›Ð?Ð¦Ð?Ð¯ ÑÑ‚Ð°Ñ€Ñ‚Ð¾Ð²Ñ‹Ñ… Ð·Ð½Ð°Ñ‡ÐµÐ½Ð¸Ð¹ Ð¿Ð°Ñ€Ð°Ð¼ÐµÑ‚Ñ€Ð¾Ð² Ñ‚Ñ€ÐµÐ½Ð´Ð°
///TODO Ð¿Ñ€Ð¸Ð´ÑƒÐ¼Ð°Ñ‚ÑŒ ÐºÐ°Ðº ÑÐ´ÐµÐ»Ð°Ñ‚ÑŒ Ð¾Ñ‡Ð¸Ð²Ð¸Ð´Ð½ÐµÐµ.
TrendParameters currentTrendParametrs = InitStartTrendParameters_ForVelocity(
        AlphaTypeConst);

//ÐšÐ¾Ð»Ð¸Ñ‡ÐµÑÑ‚Ð²Ð¾ ÑÐ±Ñ€Ð¾ÑÐ¾Ð² ÑÐºÐ¾Ñ€Ð¾ÑÑ‚Ð¸ Ð² Ð½Ð¾Ð»ÑŒ Ð·Ð° Ð¿Ð¾ÑÐ»ÐµÐ´Ð½ÐµÐµ Ð²Ñ€ÐµÐ¼Ñ
int kLoseZero = 0;
//ÐšÐ¾Ð»Ð¸Ñ‡ÐµÑÑ‚Ð²Ð¾ Ð¾ÑˆÐ¸Ð±Ð¾Ñ‡Ð½Ñ‹Ñ… Ð¾Ð¿Ñ€ÐµÐ´ÐµÐ»ÐµÐ½Ð¸Ð¹ Ð°Ð±ÑÐ¾Ð»ÑŽÑ‚Ð½Ð¾Ð³Ð¾ Ð·Ð½Ð°Ñ‡ÐµÐ½Ð¸Ñ ÑÐºÐ¾Ñ€Ð¾ÑÑ‚Ð¸ Ð·Ð° Ð¿Ð¾ÑÐ»ÐµÐ´Ð½ÐµÐµ Ð²Ñ€ÐµÐ¼Ñ
int kLoseOther = 0;

int IndexFromBigVelocity = 0;
int IndexFromZeroVelocity = 0;

int ErrorZero[ErrorZeroArraySize] = { 0 };
int ErrorOther[ErrorOtherArraySize] = { 0 };
float ShiftOther[ShiftOtherArraySize] = { 0 };

TrendParameters CalculateTrendParameters(
        float Voptic, float Vacc, float V_middle_bigFFT, float V_3_4_bigFFT,
        float V_end_bigFFT, float Health) {

    int CorrectWhileLose = 0;

    bool CHC2 = false;
    bool CHC3 = false;
    bool CHC5 = false;
    bool CHC7 = false;

    float CurrentRangeForVal = 10; // general range for Speed Error before decision about lose
    float CurrentRangeForZero = 5; // maximum output speed when Vopt=0 is valid

    float minSpeedForOptical = 0.3;

    float maxShiftOther = 0.7f; //if speed is trying to slip from the range, alpha will be increased
    float THRESH_FOR_ErZero = 0.8f; // Threshold for mean of ErZero
    float THRESH_FOR_ErOther = 0.9f; // Threshold for mean of ErOther

    int kLoseZeroMax = NumberOfOpticCalculationsInSecond * 10;
    int kLoseOtherMax = NumberOfOpticCalculationsInSecond * 2;
    int kLoseOtherMaxForLowSpeed = NumberOfOpticCalculationsInSecond * 6;

    bool goodVoptic = false;
    if (CHC2) {
        //TODO Ð”Ð¾Ð±Ð°Ð²Ð¸Ñ‚ÑŒ Ð²ÐµÑ‚ÐºÑƒ
        //goodVoptic =

    } else {

        goodVoptic = ((Voptic != 0)
                && (((std::abs(Voptic - V_middle_bigFFT) < CurrentRangeForVal)
                        && (std::abs(V_middle_bigFFT) < 50)
                        && (std::abs(V_middle_bigFFT) > minSpeedForOptical))
                        || ((std::abs(Voptic - V_middle_bigFFT)
                                < (PartOfAllowed * std::abs(V_middle_bigFFT)))
                                && (std::abs(V_middle_bigFFT) > 50))))
                || ((Voptic == 0)
                        && (std::abs(V_middle_bigFFT) < CurrentRangeForZero));
    }

    if (goodVoptic) {
        if (CHC7) {

        } else {
            // Ð¡Ð¾Ñ…Ñ€Ð°Ð½ÑÐµÐ¼ Ð¾ÑˆÐ¸Ð±ÐºÑƒ Ð² ShiftOther
            addToTheEndOfArrayAndShift(
                    ShiftOther, ShiftOtherArraySize,
                    (Voptic - V_middle_bigFFT));
        }

        float NewCurrentRangeForZero;
        if (CHC3) {
            NewCurrentRangeForZero = 2.5f;
        } else {
            NewCurrentRangeForZero = CurrentRangeForZero;
        }

        float ShiftOtherMean = math::mean(ShiftOther, ShiftOtherArraySize);
        bool errorIsBig = ((std::abs(ShiftOtherMean) > maxShiftOther)
                && (std::abs(V_middle_bigFFT) > NewCurrentRangeForZero));

        float AlphaAdd = 0;
        if (errorIsBig) {
            AlphaAdd = 0.07;
        } else {
            AlphaAdd = 0;
        }

        if (std::abs(V_middle_bigFFT) < CurrentRangeForZero) {
            currentTrendParametrs.alpha = (0.07f + kLoseZero / 1000.0)
                    + AlphaAdd;
        } else {
            currentTrendParametrs.alpha = (AlphaTypeConst + kLoseZero / 1000.0)
                    + AlphaAdd;
        }
        currentTrendParametrs.alpha1 = LinearTrend::calculateAlpha1(
                currentTrendParametrs.alpha);

        addToTheEndOfArrayAndShift(ErrorZero, ErrorZeroArraySize, 0);
        addToTheEndOfArrayAndShift(ErrorOther, ErrorOtherArraySize, 0);
        kLoseZero = 0;
        kLoseOther = 0;
    } else if ((Voptic == 0) || ((Health < 0.5f) && CHC5)) {
        /// ÐžÑˆÐ¸Ð±ÐºÐ° ÑÐ±Ñ€Ð¾Ñ Ð² Ð½Ð¾Ð»ÑŒ
        /// TODO Ð¾Ñ‚Ñ€ÐµÑ„Ð°ÐºÑ‚Ð¾Ñ€Ð¸Ñ‚ÑŒ
        currentTrendParametrs.alpha = 0;
        currentTrendParametrs.alpha1 = 0;
        addToTheEndOfArrayAndShift(ErrorZero, ErrorZeroArraySize, 1);
        addToTheEndOfArrayAndShift(ErrorOther, ErrorOtherArraySize, 0);
        kLoseZero++;

        float ErrorZeroMean = math::mean(ErrorZero, ErrorZeroArraySize);

        //TODO Ð´Ð¾Ð±Ð°Ð²Ð¸Ñ‚ÑŒ Ð¾Ð±Ð¾Ð·Ð½Ð°Ñ‡ÐµÐ½Ð¸Ðµ Ð»Ð¾Ð³Ð¸Ñ‡ÐµÑÐºÐ¾Ð³Ð¾ Ð²Ñ‹Ñ€Ð°Ð¶ÐµÐ½Ð¸Ñ
        if ((std::abs(ErrorZeroMean) > THRESH_FOR_ErZero)
                || (kLoseZero > kLoseZeroMax)) {

            reset(ErrorZero, ErrorZeroArraySize);
            kLoseZero = 0;
            CorrectWhileLose = 1;
        }
    } else {
        //TODO Ð´Ð¾Ð±Ð°Ð²Ð¸Ñ‚ÑŒ Ð¾Ð±Ð¾Ð·Ð½Ð°Ñ‡ÐµÐ½Ð¸Ðµ ÐºÐµÐ¹ÑÐ°
        currentTrendParametrs.alpha = 0;
        currentTrendParametrs.alpha1 = 0;
        addToTheEndOfArrayAndShift(ErrorZero, ErrorZeroArraySize, 0);

        if (Voptic > V_middle_bigFFT) {
            addToTheEndOfArrayAndShift(ErrorOther, ErrorOtherArraySize, 1);
        } else {
            addToTheEndOfArrayAndShift(ErrorOther, ErrorOtherArraySize, -1);
        }
        kLoseOther++;

        float ErrorOtherMean = math::mean(ErrorOther, ErrorOtherArraySize);
        if ((std::abs(ErrorOtherMean) > THRESH_FOR_ErOther)
                || ((kLoseOther > kLoseOtherMax) && (std::abs(V_middle_bigFFT) > minSpeedForOptical))
                || ((kLoseOther > kLoseOtherMaxForLowSpeed) && (std::abs(V_middle_bigFFT) < minSpeedForOptical))) {

            for (int i = 0; i < ErrorOtherArraySize; i++) {
                ErrorOther[i] = 0;
            }
            kLoseOther = 0;
            CorrectWhileLose = 2;
        }
    }
    int MaxIndexFromBigVelocity = NumberOfOpticCalculationsInSecond * 6;
    int MaxIndexFromZeroVelocity = NumberOfOpticCalculationsInSecond * 6;

    float LowSpeed = 2.0f;
    float alphaZero = 0.25f;

    //ÐšÐ¾Ñ€Ñ€ÐµÐºÑ†Ð¸Ñ alpha Ð´Ð»Ñ Ñ‚Ð¾Ñ€Ð¼Ð¾Ð·Ð½Ð¾Ð³Ð¾ Ñ‚ÐµÑÑ‚Ð°
    if ((std::abs(V_middle_bigFFT) < LowSpeed)
            && (IndexFromBigVelocity < MaxIndexFromBigVelocity)
            || (IndexFromZeroVelocity < MaxIndexFromZeroVelocity)) {
        if (Voptic == 0) {
            currentTrendParametrs.alpha = alphaZero;
            currentTrendParametrs.alpha1 = LinearTrend::calculateAlpha1(
                    AlphaTypeConst);
        }
    }

    // ÐšÐ¾Ñ€Ñ€ÐµÐºÑ†Ð¸Ñ (ÑÐ±Ñ€Ð¾Ñ)
    if (CorrectWhileLose > 0) {
        currentTrendParametrs.alpha = 0;
        currentTrendParametrs.alpha1 = 0;
        reset(ErrorZero, ErrorZeroArraySize);
        reset(ErrorOther, ErrorOtherArraySize);
    }

    IndexFromBigVelocity++;
    IndexFromZeroVelocity++;

    if (IndexFromZeroVelocity > 10000) {
        IndexFromZeroVelocity = 2 * MaxIndexFromZeroVelocity;
    }

    if (IndexFromBigVelocity > 10000) {
        IndexFromBigVelocity = 2 * MaxIndexFromBigVelocity;
    }

    if ((IndexFromBigVelocity < MaxIndexFromBigVelocity)
            && (V_3_4_bigFFT * V_end_bigFFT) < 0) {
        IndexFromZeroVelocity = 0;
    }

    if ((std::abs(V_end_bigFFT) < LowSpeed)
            && (std::abs(V_3_4_bigFFT) > LowSpeed)) {
        IndexFromBigVelocity = 0;
    }

    currentTrendParametrs.correctWhileLoose = CorrectWhileLose;

    return currentTrendParametrs;

}

/**
 * ÐžÐ±Ð½ÑƒÐ»ÑÐµÑ‚ Ð¼Ð°ÑÑÐ¸Ð²
 * @param array ÑƒÐºÐ°Ð·Ð°Ñ‚ÐµÐ»ÑŒ Ð½Ð° Ð¼Ð°ÑÑÐ¸Ð²
 * @param size Ñ€Ð°Ð·Ð¼ÐµÑ€
 */
void reset(int array[], int size) {
    std::fill(array, array + size, 0);
}

void addToTheEndOfArrayAndShift(int *array, int size, int valueToAdd) {
    std::rotate(&array[0], &array[1], &array[size]);  // ÑÐ´Ð²Ð¸Ð³ Ð²Ð»ÐµÐ²Ð¾
    array[size - 1] = valueToAdd;
}

void addToTheEndOfArrayAndShift(float *array, int size, float valueToAdd) {
    std::rotate(&array[0], &array[1], &array[size]);  // ÑÐ´Ð²Ð¸Ð³ Ð²Ð»ÐµÐ²Ð¾
    array[size - 1] = valueToAdd;
}