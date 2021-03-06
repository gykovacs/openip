#ifndef GENERALDECISIONTREECLASSIFIER_H
#define GENERALDECISIONTREECLASSIFIER_H
#include <openipDS/DataSet.h>
#include <openipDS/mathFunctions.h>

#include <openipML/BinaryDecisionTreeNode.h>
#include <openipML/Classifier.h>
#include <openipML/DecisionTreeClassifier.h>

namespace openip
{
    template<typename T>
    class GeneralDecisionTreeClassifier: public Classifier<T>
    {
    public:
        using Classifier<T>::probabilities;
        using Classifier<T>::saveModelToFile;
        using Classifier<T>::openModelFromFile;

        /**
         * default constructor
         */
        GeneralDecisionTreeClassifier();

        /**
         * copy constructor
         * @param c instance to copy
         */
        GeneralDecisionTreeClassifier(const GeneralDecisionTreeClassifier& c);

        /**
         * destructor
         */
        ~GeneralDecisionTreeClassifier();

        /**
         * training method of the classifier
         * @param data training dataset
         * @param mask only those elements of the dataset are trained, which have a corresponding non-zero element in the mask
         * @return non-zero in case of errors, zero otherwise
         */
        virtual int train(DataSet<T>& data, Vector<unsigned char>* mask= NULL);

        /**
         * classifies the parameter vector
         * @param c vector to classify
         * @param probabilities probability estimates
         * @return predicted class index of the vector
         */
        virtual int classify(Vector<T>& c, Vector<float>* probabilities= NULL);

        void determineSplit(DataSet<T>& data, Vector<unsigned char>& mask, int& attribute, T& threshold);

        float informationGain(DataSet<T>& data, Vector<unsigned char>& mask, int coordinate, T& threshold);

        void buildTree(BinaryDecisionTreeNode<T>* &root, DataSet<T>& data, Vector<unsigned char>& mask);

        int getLabel(BinaryDecisionTreeNode<T>* node, Vector<T>& c);

        void deleteTree(BinaryDecisionTreeNode<T>** root);

        /**
         * saves the model to filestream
         * @param file output filestream
         * @return non-zero in case of errors, zero otherwise
         */
        virtual int saveModelToFile(std::ofstream& file);

        virtual int writeTreeIntoFile(BinaryDecisionTreeNode<T>* p, std::ofstream& file);

        /**
         * opens the model from filestream
         * @param file input filestream
         * @return non-zero in case of errors, zero otherwise
         */
        virtual int openModelFromFile(std::ifstream& file);

        virtual int readTreeFromFile(BinaryDecisionTreeNode<T>** p, std::ifstream& file);

        BinaryDecisionTreeNode<T>* root;

        T oldThreshold;
        int oldCoordinate;
    };

    template<typename T>
    GeneralDecisionTreeClassifier<T>::GeneralDecisionTreeClassifier()
    : Classifier<T>()
    {
        root= NULL;
    }

    template<typename T>
    GeneralDecisionTreeClassifier<T>::GeneralDecisionTreeClassifier(const GeneralDecisionTreeClassifier& c)
    : Classifier<T>(c)
    {
        root= c.root;
    }

    template<typename T>
    GeneralDecisionTreeClassifier<T>::~GeneralDecisionTreeClassifier()
    {
        deleteTree(&root);
    }

    template<typename T>
    void GeneralDecisionTreeClassifier<T>::deleteTree(BinaryDecisionTreeNode<T> **root)
    {
        //printf("delete tree %p\n", *root); fflush(stdout);
        if ( *root )
        {
            if ( (*root)->childrenNodes[0] )
                deleteTree(&((*root)->childrenNodes[0]));
            if ( (*root)->childrenNodes[1] )
                deleteTree(&((*root)->childrenNodes[1]));
            delete *root;
            *root= NULL;
        }
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::train(DataSet<T>& data, Vector<unsigned char>* mask)
    {
        printf("train gdtc %lu %lu %lu %lu\n", data.classNames.size(), data.featureNames.size(), data.preprocessingNames.size(), data.size()); fflush(stdout);
        deleteTree(&root);

        this->classNames= data.classNames;
        this->featureNames= data.featureNames;
        this->preprocessingNames= data.preprocessingNames;

        oldThreshold= FLT_MAX;
        oldCoordinate= -1;

        if ( mask )
        {
            DataSet<T> tmp;
            for ( unsigned int i= 0; i < mask->size(); ++i )
                if ( (*mask)(i) )
                    tmp.push_back(data(i));
            Vector<unsigned char> tmpMask(tmp.size());
            tmpMask= 1;

            tmp.init();

            buildTree(root, tmp, tmpMask);
        }
        else
        {
            Vector<unsigned char> tmpMask(data.size());
            tmpMask= 1;

            data.init();

            buildTree(root, data, tmpMask);
        }
        //printf("trining finished "); fflush(stdout);
        return 0;
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::classify(Vector<T>& c, Vector<float>* probabilities)
    {
        if ( probabilities == NULL )
            probabilities= &(this->probabilities);
        return getLabel(root, c);
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::getLabel(BinaryDecisionTreeNode<T>* node, Vector<T>& c)
    {
        if ( node->terminal )
            return node->label;
        if ( c(node->attribute) < node->threshold )
            return getLabel(node->left(), c);
        else
            return getLabel(node->right(), c);
    }

    template<typename T>
    void GeneralDecisionTreeClassifier<T>::buildTree(BinaryDecisionTreeNode<T>* &root, DataSet<T>& data, Vector<unsigned char>& mask)
    {
        int firstLabel= -1;
        Vector<T> firstValue;
        bool terminal= true;
        bool unique= true;
        for ( unsigned int i= 0; i < mask.size(); ++i )
        {
            if ( mask(i) && firstLabel == -1 )
            {
                firstLabel= data(i).classIndex;
                firstValue= data(i);
            }
            if ( firstLabel != -1 && mask(i) && data(i).classIndex != firstLabel)
            {
                terminal= false;
            }
            if ( firstLabel != -1 && mask(i) && data(i) != firstValue )
            {
                unique= false;
            }
        }

        if ( unique )
        {
            Vector<int> labels(data.numberOfClasses);
            for ( unsigned int i= 0; i < mask.size(); ++i )
                if ( mask(i) )
                    labels(data(i).classIndex)++;
            int maxLabel= 0;
            for ( unsigned int i= 0; i < labels.size(); ++i )
                if ( labels(maxLabel) < labels(i) )
                    maxLabel= i;
            root= new BinaryDecisionTreeNode<T>(true, maxLabel);
        }
        else if ( terminal )
        {
            root= new BinaryDecisionTreeNode<T>(true, firstLabel);
        }
        else
        {
            int attribute;
            T threshold;
            determineSplit(data, mask, attribute, threshold);

            if ( attribute == oldCoordinate && threshold == oldThreshold )
            {
                Vector<int> labels(data.numberOfClasses);
                for ( unsigned int i= 0; i < mask.size(); ++i )
                    if ( mask(i) )
                        labels(data(i).classIndex)++;
                int maxLabel= 0;
                for ( unsigned int i= 0; i < labels.size(); ++i )
                    if ( labels(maxLabel) < labels(i) )
                        maxLabel= i;
                root= new BinaryDecisionTreeNode<T>(true, labels(maxLabel));

                return;
            }
            oldCoordinate= attribute;
            oldThreshold= threshold;

            root= new BinaryDecisionTreeNode<T>(false, attribute, threshold);

            Vector<unsigned char> maskLeft(data.size());
            Vector<unsigned char> maskRight(data.size());
            int left= 0, right= 0;
            for ( unsigned int i= 0; i < mask.size(); ++i )
                if ( mask(i) && data(i)(attribute) <= threshold )
                {
                    maskLeft(i)= 1;
                    maskRight(i)= 0;
                    ++left;
                }
                else if ( mask(i) && data(i)(attribute) > threshold )
                {
                    maskRight(i)= 1;
                    maskLeft(i)= 0;
                    ++right;
                }
                else
                {
                    maskLeft(i)= maskRight(i)= 0;
                }
            buildTree(root->childrenNodes[0], data, maskLeft);
            buildTree(root->childrenNodes[1], data, maskRight);
        }

    }

    template<typename T>
    void GeneralDecisionTreeClassifier<T>::determineSplit(DataSet<T>& data, Vector<unsigned char>& mask, int& attribute, T& threshold)
    {
        //printf("f"); fflush(stdout);
        int bestAttribute= 0;
        float bestThreshold= 0;
        float bestIG= FLT_MAX;

        float tmpThreshold;
        float tmpIG;

        for ( unsigned int i= 0; i < data(0).size(); ++i )
        {
            //printf("ig "); fflush(stdout);
            tmpIG= informationGain(data, mask, i, tmpThreshold);
            //printf("eig "); fflush(stdout);
            if ( tmpIG < bestIG )
            {
                bestIG= tmpIG;
                bestAttribute= i;
                bestThreshold= tmpThreshold;
            }
        }
        attribute= bestAttribute;
        threshold= bestThreshold;
        //printf("g"); fflush(stdout);
    }

    template<typename T>
    float GeneralDecisionTreeClassifier<T>::informationGain(DataSet<T>& data, Vector<unsigned char>& mask, int coordinate, T& threshold)
    {
        int numberOfVectors= 0;
        for ( unsigned int i= 0; i < mask.size(); ++i )
            if ( mask(i) )
                ++numberOfVectors;
        //printf("%d ", numberOfVectors); fflush(stdout);

        Vector<SortingRecord> sample(numberOfVectors);
        Vector<int> classFrequencies(data.numberOfClasses);
        classFrequencies= 0;
        int s= 0;
        for ( unsigned int i= 0; i < data.size(); ++i )
        {
            if ( mask(i) )
            {
                sample(s++)= SortingRecord(data(i)(coordinate), data(i).classIndex);
                classFrequencies(data(i).classIndex)++;
            }
        }

        sort(sample.begin(), sample.end(), sortByData);

        Vector<int> p;
        Vector<int> pleft;
        Vector<int> pright;

        p.resize(2);
        pleft.resize(2);
        pright.resize(2);

        p= 0;
        pleft= 0;
        pright= 0;

        for ( unsigned int i= 0; i < sample.size(); ++i )
            p(sample(i).classIndex)++;

        pright= p;
        float parentGini= 1.0 - p(0)/numberOfVectors*p(0)/numberOfVectors - p(1)/numberOfVectors*p(1)/numberOfVectors;

        float left;
        float right;
        float ginileft;
        float giniright;
        float igain;
        float bestidx= 0;
        float bestigain= FLT_MAX;
        for ( int i= 0; i < int(sample.size()); ++i )
        {
            pleft(sample(i).classIndex)++;
            pright(sample(i).classIndex)--;

            if ( i != int(sample.size())-1 && sample(i).data == sample(i+1).data )
                 continue;

            left= pleft(0) + pleft(1);
            right= pright(0) + pright(1);

            ginileft= 1.0 - pleft(0)/left*pleft(0)/left - pleft(1)/left*pleft(1)/left;
            giniright= 1.0 - pright(0)/right*pright(0)/left - pright(1)/right*pright(1)/right;

            igain= parentGini - left/numberOfVectors*ginileft - right/numberOfVectors*giniright;

            if ( igain < bestigain )
            {
                bestigain= igain;
                bestidx= i;
            }
        }

        threshold= sample(bestidx).data;
        return bestigain;
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::saveModelToFile(std::ofstream& file)
    {
        file << "GeneralDecisionTreeClassifier" << std::endl;
        this->Classifier<T>::saveModelToFile(file);

        writeTreeIntoFile(root, file);

        return 0;
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::writeTreeIntoFile(BinaryDecisionTreeNode<T>* p, std::ofstream& file)
    {
        if ( p == 0 )
        {
            file << "0" << std::endl;
        }
        else
        {
            file << "1 " << p->attribute << " " << p->label << " " << p->numberOfChildren << " " << p->terminal << " " << p->threshold << std::endl;
            writeTreeIntoFile(p->left(), file);
            writeTreeIntoFile(p->right(), file);
        }

        return 0;
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::openModelFromFile(std::ifstream& file)
    {
        std::string tmp;

        tprintf("dtc: open model from file\n");
        getline(file, tmp);

        if ( tmp.compare("GeneralDecisionTreeClassifier") == 0 )
        {
            tprintf("dtc: open: in if\n");
            this->Classifier<T>::openModelFromFile(file);
            readTreeFromFile(&root, file);

            return 0;
        }

        probabilities.resize(this->classNames.size());
        tprintf("decision tree read successfully\n");
        return 1;
    }

    template<typename T>
    int GeneralDecisionTreeClassifier<T>::readTreeFromFile(BinaryDecisionTreeNode<T>** p, std::ifstream &file)
    {
        if ( file.eof() )
        {
            *p= NULL;
            return 0;
        }

        int token;
        file >> token;

        int attribute;
        int label;
        bool terminal;
        int numberOfChildren;
        T threshold;

        if ( token == 0 )
        {
            *p= 0;
            return 0;
        }
        else
        {
            file >> attribute;
            file >> label;
            file >> numberOfChildren;
            file >> terminal;
            file >> threshold;

            *p= new BinaryDecisionTreeNode<T>(terminal, attribute, threshold);
            (*p)->label= label;
            (*p)->numberOfChildren= numberOfChildren;

            readTreeFromFile(&((*p)->childrenNodes[0]), file);
            readTreeFromFile(&((*p)->childrenNodes[1]), file);
        }

        return 0;
    }

}

#endif // GENERALDECISIONTREECLASSIFIER_H
