#include <stdio.h>
#include <vector>
#include <array>
#include <random>
#include <string>

#define VERBOSE() false  // if true, shows more comparisons
#define DETERMINISTIC() true
static const int c_numItems = 10;
static const int c_numRollsShow = 80;
static const int c_numRollsHistogram[] = { 10, 100, 1000, 10000, 100000, 1000000 };
static const char c_baseCharacter = '0';

static const int c_numRollsTotal = c_numRollsHistogram[sizeof(c_numRollsHistogram) / sizeof(c_numRollsHistogram[0]) - 1];
static const float c_goldenRatioConjugate = 0.61803398875f;
static const float c_piFract = 0.14159265359f;
static const float c_sqrt2Fract = 0.41421356237f;

typedef std::vector<std::vector<std::string>> CSV;
typedef std::array<float, 2> Vec2;

std::mt19937 GetRNG()
{
#if DETERMINISTIC()
    std::mt19937 rng;
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif
    return rng;
}

float fract(float f)
{
    return f - floor(f);
}

int min(int a, int b)
{
    return a <= b ? a : b;
}

int FloatToItem(float f, int numItems)
{
    // if f is in [0,1], remaps to [0, numItems-1]
    return min(int(f * float(numItems)), numItems - 1);
}

int FloatToWeightedItem(float f, const std::vector<float>& normalizedItemWeights)
{
    for (size_t i = 0; i < normalizedItemWeights.size(); ++i)
    {
        f -= normalizedItemWeights[i];
        if (f <= 0.0f)
            return int(i);
    }

    return int(normalizedItemWeights.size() - 1);
}

void ShowSequence(const char* label, const std::vector<int>& sequence, int count)
{
    printf("  %s:\n    ", label);

    for (int i = 0; i < count; ++i)
        printf("%c", c_baseCharacter + sequence[i]);

    printf("\n\n");
}

void AddHistogram(CSV& csv, const char* label, const std::vector<int>& sequence, int count)
{
    std::vector<float> histogram;
    histogram.resize(c_numItems, 0.0f);

    for (int i = 0; i < count; ++i)
        histogram[sequence[i]] += 1.0f;

    for (float& f : histogram)
        f /= float(count);

    size_t rowIndex = csv.size();
    csv.resize(rowIndex + 1);
    std::vector<std::string>& row = csv[rowIndex];

    row.push_back(label);

    char buffer[1024];
    for (int i = 0; i < c_numItems; ++i)
    {
        sprintf_s(buffer, "%f", histogram[i]);
        row.push_back(buffer);
    }
}

void SaveCSV(const CSV& csv, const char* type, int count)
{
    char fileName[256];
    sprintf_s(fileName, "out/histogram_%s_%i.csv", type, count);

    FILE* file = nullptr;
    fopen_s(&file, fileName, "w+t");

    for (auto& row : csv)
    {
        bool firstCell = true;
        for (auto& cell : row)
        {
            fprintf(file, "%s\"%s\"", firstCell ? "" : ",", cell.c_str());
            firstCell = false;
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

Vec2 R2Additive(const Vec2& input)
{
    // generalized golden ratio, from:
    // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    const float g = 1.32471795724474602596f;
    const float a1 = 1.0f / g;
    const float a2 = 1.0f / (g * g);

    return Vec2{
        fract(input[0] + a1),
        fract(input[1] + a2),
    };
}

Vec2 R2(int index)
{
    // generalized golden ratio, from:
    // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    const float g = 1.32471795724474602596f;
    const float a1 = 1.0f / g;
    const float a2 = 1.0f / (g * g);
    return Vec2{fract(a1 * float(index)), fract(a2 * float(index))};
}

static size_t Ruler(size_t n)
{
    size_t ret = 0;
    while (n != 0 && (n & 1) == 0)
    {
        n /= 2;
        ++ret;
    }
    return ret;
}

void Sobol(std::vector<Vec2>& values, size_t numValues)
{
    // x axis
    values.resize(numValues);
    size_t sampleInt = 0;
    for (size_t i = 0; i < numValues; ++i)
    {
        size_t ruler = Ruler(i + 1);
        size_t direction = size_t(size_t(1) << size_t(31 - ruler));
        sampleInt = sampleInt ^ direction;
        values[i][0] = float(sampleInt) / std::pow(2.0f, 32.0f);
    }

    // y axis
    // Code adapted from http://web.maths.unsw.edu.au/~fkuo/sobol/
    // uses numbers: new-joe-kuo-6.21201

    // Direction numbers
    std::vector<size_t> V;
    V.resize((size_t)ceil(log((double)numValues + 1) / log(2.0)));  //+1 because we are skipping index 0
    V[0] = size_t(1) << size_t(31);
    for (size_t i = 1; i < V.size(); ++i)
        V[i] = V[i - 1] ^ (V[i - 1] >> 1);

    // Samples
    sampleInt = 0;
    for (size_t i = 0; i < numValues; ++i) {
        size_t ruler = Ruler(i + 1);
        sampleInt = sampleInt ^ V[ruler];
        values[i][1] = float(sampleInt) / std::pow(2.0f, 32.0f);
    }
}

struct AliasTableEntry
{
    float probability;
    int aliasIndex;
};

typedef std::vector<AliasTableEntry> AliasTable;

int SampleAliasTable(const AliasTable& aliasTable, float x, float y)
{
    int column = FloatToItem(x, int(aliasTable.size()));
    const AliasTableEntry& entry = aliasTable[column];
    if (entry.probability >= 1.0f || y <= entry.probability)
        return column;
    else
        return entry.aliasIndex;
}

AliasTable MakeAliasTable(const std::vector<float>& itemWeights)
{
    // uses the stable vose alias method from https://www.keithschwarz.com/darts-dice-coins/
    AliasTable aliasTable(itemWeights.size());

    // make the small and large item lists
    struct Item
    {
        int index;
        float probability;
    };
    std::vector<Item> largeItems, smallItems;
    for (size_t index = 0; index < itemWeights.size(); ++index)
    {
        Item newItem;
        newItem.index = int(index);
        newItem.probability = itemWeights[index] * float(itemWeights.size());

        if (newItem.probability < 1.0f)
            smallItems.push_back(newItem);
        else
            largeItems.push_back(newItem);
    }

    while (!largeItems.empty() && !smallItems.empty())
    {
        Item smallItem = smallItems.back();
        smallItems.pop_back();

        Item largeItem = largeItems.back();
        largeItems.pop_back();

        aliasTable[smallItem.index].probability = smallItem.probability;
        aliasTable[smallItem.index].aliasIndex = largeItem.index;

        largeItem.probability = (largeItem.probability + smallItem.probability) - 1.0f;
        if (largeItem.probability < 1.0f)
            smallItems.push_back(largeItem);
        else
            largeItems.push_back(largeItem);
    }

    for (const Item& item : largeItems)
    {
        aliasTable[item.index].probability = 1.0f;
        aliasTable[item.index].aliasIndex = -1;
    }

    for (const Item& item : smallItems)
    {
        aliasTable[item.index].probability = 1.0f;
        aliasTable[item.index].aliasIndex = -1;
    }

    return aliasTable;
}

int main(int argc, char** argv)
{
    // do non weighted tests
    {
        std::mt19937 rng = GetRNG();
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        std::vector<int> seq_Sequential;
        std::vector<int> seq_WhiteNoise;
        std::vector<int> seq_GoldenRatio;
        std::vector<int> seq_OneMinusGoldenRatio;
        std::vector<int> seq_Pi;
        std::vector<int> seq_OneMinusPi;
        std::vector<int> seq_sqrt2;

        // Note: these could start at any value between 0 and 1.
        float goldenRatioValue = 0.0f;
        float oneMinusGoldenRatioValue = 0.0f;
        float piValue = 0.0f;
        float oneMinusPiValue = 0.0f;
        float sqrt2Value = 0.0f;

        float sequentialDelta = 1.0f / float(c_numItems);
        float sequentialValue = sequentialDelta / 2.0f; // to help numerical problems from it being right on the edge.

        // generate the items, using each sequence type
        for (int i = 0; i < c_numRollsTotal; ++i)
        {
            sequentialValue = fract(sequentialValue + sequentialDelta);
            seq_Sequential.push_back(FloatToItem(sequentialValue, c_numItems));

            float whiteNoiseValue = dist(rng);
            seq_WhiteNoise.push_back(FloatToItem(whiteNoiseValue, c_numItems));

            goldenRatioValue = fract(goldenRatioValue + c_goldenRatioConjugate);
            seq_GoldenRatio.push_back(FloatToItem(goldenRatioValue, c_numItems));

            oneMinusGoldenRatioValue = fract(oneMinusGoldenRatioValue + (1.0f - c_goldenRatioConjugate));
            seq_OneMinusGoldenRatio.push_back(FloatToItem(oneMinusGoldenRatioValue, c_numItems));

            piValue = fract(piValue + c_piFract);
            seq_Pi.push_back(FloatToItem(piValue, c_numItems));

            oneMinusPiValue = fract(oneMinusPiValue + (1.0f - c_piFract));
            seq_OneMinusPi.push_back(FloatToItem(oneMinusPiValue, c_numItems));

            sqrt2Value = fract(sqrt2Value + c_sqrt2Fract);
            seq_sqrt2.push_back(FloatToItem(sqrt2Value, c_numItems));
        }

        // show items 
        printf("=================== Unweighted ===================\n\n");
        ShowSequence("Sequential", seq_Sequential, c_numRollsShow);
        ShowSequence("White Noise", seq_WhiteNoise, c_numRollsShow);
        ShowSequence("Golden Ratio", seq_GoldenRatio, c_numRollsShow);
        ShowSequence("Pi", seq_Pi, c_numRollsShow);
        ShowSequence("Sqrt2", seq_sqrt2, c_numRollsShow);
#if VERBOSE()
        ShowSequence("One Minus Golden Ratio", seq_OneMinusGoldenRatio, c_numRollsShow);
        ShowSequence("One Minus Pi", seq_OneMinusPi, c_numRollsShow);
#endif

        // write out histograms to csvs at each step of c_numRollsHistogram
        for (size_t i = 0; i < sizeof(c_numRollsHistogram) / sizeof(c_numRollsHistogram[0]); ++i)
        {
            CSV csv;
            AddHistogram(csv, "Sequential", seq_Sequential, c_numRollsHistogram[i]);
            AddHistogram(csv, "White Noise", seq_WhiteNoise, c_numRollsHistogram[i]);
            AddHistogram(csv, "Golden Ratio", seq_GoldenRatio, c_numRollsHistogram[i]);
            AddHistogram(csv, "Pi", seq_Pi, c_numRollsHistogram[i]);
            AddHistogram(csv, "Sqrt2", seq_sqrt2, c_numRollsHistogram[i]);
#if VERBOSE()
            AddHistogram(csv, "One Minus Golden Ratio", seq_OneMinusGoldenRatio, c_numRollsHistogram[i]);
            AddHistogram(csv, "One Minus Pi", seq_OneMinusPi, c_numRollsHistogram[i]);
#endif
            SaveCSV(csv, "unweighted", c_numRollsHistogram[i]);
        }
    }

    // do weighted tests
    {
        // calculate some weights for items and then normalize them so they add up to 1 and are a pmf
        std::vector<float> itemWeights;
        float weightTotal = 0.0f;
        for (int i = 1; i <= c_numItems; ++i)
        {
            itemWeights.push_back(float(i));
            weightTotal += float(i);
        }
        for (float& f : itemWeights)
            f /= weightTotal;
        float smallestWeight = itemWeights[0];

        AliasTable aliasTable = MakeAliasTable(itemWeights);

        std::mt19937 rng = GetRNG();
        std::mt19937 rng2 = GetRNG();
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        std::vector<int> seq_Sequential;
        std::vector<int> seq_WhiteNoise;
        std::vector<int> seq_GoldenRatio;
        std::vector<int> seq_OneMinusGoldenRatio;
        std::vector<int> seq_Pi;
        std::vector<int> seq_OneMinusPi;
        std::vector<int> seq_sqrt2;

        std::vector<int> seq_AliasWhiteNoise;
        std::vector<int> seq_AliasR2;
        std::vector<int> seq_AliasR2Additive;
        std::vector<int> seq_AliasGrSqrt2;
        std::vector<int> seq_AliasSobol;

        // Generate sobol samples
        std::vector<Vec2> sobolSamples;
        Sobol(sobolSamples, c_numRollsTotal);

        // Note: these could start at any value between 0 and 1.
        float goldenRatioValue = 0.0f;
        float oneMinusGoldenRatioValue = 0.0f;
        float piValue = 0.0f;
        float oneMinusPiValue = 0.0f;
        float sqrt2Value = 0.0f;
        Vec2 R2Value = { 0.0f, 0.0f };

        float sequentialDelta = smallestWeight;
        float sequentialValue = sequentialDelta / 2.0f; // to help numerical problems from it being right on the edge.

        // generate the items, using each sequence type
        for (int i = 0; i < c_numRollsTotal; ++i)
        {
            sequentialValue = fract(sequentialValue + sequentialDelta);
            seq_Sequential.push_back(FloatToWeightedItem(sequentialValue, itemWeights));

            float whiteNoiseValue = dist(rng);
            seq_WhiteNoise.push_back(FloatToWeightedItem(whiteNoiseValue, itemWeights));

            goldenRatioValue = fract(goldenRatioValue + c_goldenRatioConjugate);
            seq_GoldenRatio.push_back(FloatToWeightedItem(goldenRatioValue, itemWeights));

            oneMinusGoldenRatioValue = fract(oneMinusGoldenRatioValue + (1.0f - c_goldenRatioConjugate));
            seq_OneMinusGoldenRatio.push_back(FloatToWeightedItem(oneMinusGoldenRatioValue, itemWeights));

            piValue = fract(piValue + c_piFract);
            seq_Pi.push_back(FloatToWeightedItem(piValue, itemWeights));

            oneMinusPiValue = fract(oneMinusPiValue + (1.0f - c_piFract));
            seq_OneMinusPi.push_back(FloatToWeightedItem(oneMinusPiValue, itemWeights));

            sqrt2Value = fract(sqrt2Value + c_sqrt2Fract);
            seq_sqrt2.push_back(FloatToWeightedItem(sqrt2Value, itemWeights));

            // Alias table methods begin

            seq_AliasWhiteNoise.push_back(SampleAliasTable(aliasTable, dist(rng2), dist(rng2)));

            Vec2 R2Value = R2(i);
            seq_AliasR2.push_back(SampleAliasTable(aliasTable, R2Value[0], R2Value[1]));

            R2Value = R2Additive(R2Value);
            seq_AliasR2Additive.push_back(SampleAliasTable(aliasTable, R2Value[0], R2Value[1]));

            seq_AliasGrSqrt2.push_back(SampleAliasTable(aliasTable, goldenRatioValue, sqrt2Value));

            seq_AliasSobol.push_back(SampleAliasTable(aliasTable, sobolSamples[i][0], sobolSamples[i][1]));
        }

        // show items 
        printf("=================== Weighted ===================\n\n");
        ShowSequence("Sequential", seq_Sequential, c_numRollsShow);
        ShowSequence("White Noise", seq_WhiteNoise, c_numRollsShow);
        ShowSequence("Golden Ratio", seq_GoldenRatio, c_numRollsShow);
        ShowSequence("Pi", seq_Pi, c_numRollsShow);
        ShowSequence("Sqrt2", seq_sqrt2, c_numRollsShow);
        ShowSequence("Alias White Noise", seq_AliasWhiteNoise, c_numRollsShow);
        ShowSequence("Alias R2", seq_AliasR2, c_numRollsShow);
        ShowSequence("Alias R2 (Additive)", seq_AliasR2Additive, c_numRollsShow);
        ShowSequence("Alias GR / Sqrt2", seq_AliasGrSqrt2, c_numRollsShow);
        ShowSequence("Alias Sobol", seq_AliasSobol, c_numRollsShow);
#if VERBOSE()
        ShowSequence("One Minus Golden Ratio", seq_OneMinusGoldenRatio, c_numRollsShow);
        ShowSequence("One Minus Pi", seq_OneMinusPi, c_numRollsShow);
#endif

        // write out histograms to csvs at each step of c_numRollsHistogram
        for (size_t i = 0; i < sizeof(c_numRollsHistogram) / sizeof(c_numRollsHistogram[0]); ++i)
        {
            CSV csv;

        // TODO: maybe have a second group in the histogram for error ?

            // show the target value.
            {
                size_t rowIndex = csv.size();
                csv.resize(rowIndex + 1);
                std::vector<std::string>& row = csv[rowIndex];

                row.push_back("Weights");

                char buffer[1024];
                for (float f : itemWeights)
                {
                    sprintf_s(buffer, "%f", f);
                    row.push_back(buffer);
                }
            }

            AddHistogram(csv, "Sequential", seq_Sequential, c_numRollsHistogram[i]);
            AddHistogram(csv, "White Noise", seq_WhiteNoise, c_numRollsHistogram[i]);
            AddHistogram(csv, "Golden Ratio", seq_GoldenRatio, c_numRollsHistogram[i]);
            AddHistogram(csv, "Pi", seq_Pi, c_numRollsHistogram[i]);
            AddHistogram(csv, "Sqrt2", seq_sqrt2, c_numRollsHistogram[i]);
            AddHistogram(csv, "Alias White Noise", seq_AliasWhiteNoise, c_numRollsHistogram[i]);
            AddHistogram(csv, "Alias R2", seq_AliasR2, c_numRollsHistogram[i]);
            AddHistogram(csv, "Alias R2 (Additive)", seq_AliasR2Additive, c_numRollsHistogram[i]);
            AddHistogram(csv, "Alias GR / Sqrt2", seq_AliasGrSqrt2, c_numRollsHistogram[i]);
            AddHistogram(csv, "Alias Sobol", seq_AliasSobol, c_numRollsHistogram[i]);
#if VERBOSE()
            AddHistogram(csv, "One Minus Golden Ratio", seq_OneMinusGoldenRatio, c_numRollsHistogram[i]);
            AddHistogram(csv, "One Minus Pi", seq_OneMinusPi, c_numRollsHistogram[i]);
#endif
            SaveCSV(csv, "weighted", c_numRollsHistogram[i]);
        }
    }

    system("pause");

    return 0;
}
